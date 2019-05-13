#!/usr/bin/env python

import os, sys, re, filecmp
from numpy import *
import argparse
from glob import glob


def system_run(cmd):
    print cmd
    exit = os.system(cmd)
    if ( exit != 0 ):
        print 'failed to run:'
        print cmd
        sys.exit()
 
def read_pdb(pdbfile,TER=False):
    f = open(pdbfile,'r')
    pdb = []
    res_info = []
    tot_charge = 0
    for line in f:
        record = line[:6]
        if ( record != 'ATOM  ' and record != 'HETATM' ): continue
        serial = int( line[6:11] )
        atomname = line[12:16]
        altloc = line[16]
        resname = line[17:20]
        chain = line[21]
        resnum = int( line[22:26] )
        achar = line[26]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        try:
            occ = float( line[54:60] )
        except:
            occ = 1.0
        try:
            tfactor = float( line[60:66] )
        except:
            tfactor = 1.0
        try:
            segid = line[72:76]
        except:
            segid = ''
        try:
            elsymbol = line[76:78]
        except:
            elsymbol = ''
        try: 
            charge =  line[78:80] 
        except:
            charge = '0.'
        # 0:  record
        # 1:  serial
        # 2:  atomname
        # 3:  altloc
        # 4:  resname
        # 5:  chain
        # 6:  resnum
        # 7:  achar
        # 8:  x
        # 9:  y
        # 10: z
        # 11: occ
        # 12: tfactor
        # 13: segid
        # 14: elsymbol
        # 15: charge
        try:
            fix = line[85:87]
            if int(fix) == -1:
                res_info.append([atomname, chain, resnum])
        except:
            fix = " 0"
        if '+' in charge:
            tot_charge += int(charge[0])    
        elif '-' in charge:
            tot_charge -= int(charge[0])
        else:
            tot_charge += 0
            
        pdb.append( [ record, serial, atomname, altloc, resname, chain, resnum, achar, x, y, z, occ, tfactor, segid, elsymbol, charge, fix ] )
        #     1          0       1       2       3       4           5   6       7      8  9  10 11   12         13      14      15      16
    f.close()
    return pdb, res_info, tot_charge

def write_pdb(filename,pdb,renum_atom=True,hydrogen=True,renum_res=False):
    if ( isinstance(filename,file) ):
        f = filename
        file_opened = False
    else:
        f = open(filename,'w')
        file_opened = True
    serial = 1
    serial_res = 1
    prev_res = pdb[0][6]
    for p in pdb:
        t = p[:] # copy
        if p[0] == 'TER':
            f.write('TER\n')
            continue
        if ( hydrogen == False and p[2].strip()[0] == 'H' ): continue
        if ( renum_atom ):
            t = [t[0],serial] + t[2:17]
            serial += 1
        else:
            t = t[:17]
        if ( renum_res ):
            if ( prev_res != p[6] ):
                serial_res += 1
                prev_res = p[6]
            t = t[0:6]+[serial_res]+t[7:17]
        f.write('%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%3s%6s\n'%tuple(t))
    if ( file_opened ):
        f.close()

### copy h-added pdb xyz and other information into tmppdb ###
def pdb_after_addh(tmppdb,newpdb):
    tmp_pdb, res_info, tot_charge_t = read_pdb(tmppdb)
    tmp_xyz = []
    for i in tmp_pdb:
        tmp_xyz.append([i[8],i[9],i[10]])
    new_pdb, binfo, tot_charge = read_pdb(newpdb)     #can be just xyz files from cerius or pymol
    pic_atom = []
    for line in new_pdb:
        if [line[8],line[9],line[10]] not in tmp_xyz:
            line[15] = '0 '
            pic_atom.append([line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],' H',line[15],line[16]])
        else:
            if '+' in line[15] or '-' in line[15]:
                charge = line[15]
            else:
                charge = '0 '
            idx = tmp_xyz.index([line[8],line[9],line[10]])
            line = tmp_pdb[idx]
            line[15] = charge
            if [line[2],line[5],line[6]] in res_info:
                pic_atom.append([line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],'-1'])
            else:
                pic_atom.append([line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],line[16]])
    return pic_atom, tot_charge #, xyz, atom, hold

def write_input(inp_name,inp_temp,charge,multiplicity,pic_atom,tot_charge):
    ### inp_name default can be 1.inp, but first is name.input such as 9.input
    ### inp_type list which includes [small/large, level_theory, basis, opt, freq,  
    ### input_template line0: size
    ### input_template line1: level
    ### input_template line2: opt + calcfc/readfc or opt + modred + info
    ### input_template line3: freq
    ### input_template line4: guess=read
    ### input_template line5: geom=checkpoint
    ### input_template line6: iop
    ### input_template line7: scf_info
    ### input_template line8: scrf_info

    with open(inp_temp) as f:
        lines = f.readlines()

    inp = open('%s'%inp_name,'w')
    inp.write("%chk=1.chk\n")           #write check file into 1.chk
    inp.write("%nprocshared=10\n")      #use 10 processors

    if 'small' in lines[0]:
        inp.write("%mem=20GB\n")        #use memory if it is large job use 80GB
    else:
        inp.write("%mem=80GB\n")
    
    inp.write("#P %s "%lines[1].strip())      
    if bool(lines[2].strip()):
        optl = lines[2].split()
    inp.write("%s "%optl[0])
    if 'modred' in lines[2]:
        f_atom = []
        modred_info, modred_code = optl[1:3]
        pairs = modred_info.split(';')
        lmod = []
        for pair in pairs:
            ar = pair.split(',')
            lmod.append(len(ar)/2)
            for i in range(0,len(ar),2):
                for atom in pic_atom:
                    if ar[i] in atom[2] and ar[i+1] in atom[4]:
                        f_atom.append(pic_atom.index(atom)+1)
    for l in range(3,9):
        if bool(lines[l].strip()):
            inp.write("%s "%lines[l].strip())

    inp.write("\n\n")
    inp.write("info_line\n")
    inp.write("\n")
    inp.write("%d %d\n"%(charge+tot_charge,multiplicity))

    if 'check' not in lines[5]:
        ### pm7 with opt only will be relax h step/ pm7 with opt(modred) otherwise
        if 'pm7' in lines[1] and lines[2].strip() == 'opt':
            for atom in pic_atom:
                if atom[14].strip() == 'H':
                    inp.write("%4s %6s         %8.3f %8.3f %8.3f\n"%(atom[14].strip(),'0',atom[8],atom[9],atom[10])) 
                else:
                    inp.write("%4s %6s         %8.3f %8.3f %8.3f\n"%(atom[14].strip(),'-1',atom[8],atom[9],atom[10])) 
        else:
            for atom in pic_atom:
                inp.write("%4s %6s         %8.3f %8.3f %8.3f\n"%(atom[14].strip(),atom[16],atom[8],atom[9],atom[10])) 
    inp.write("\n")

    count = 0
    if 'modred' in lines[2]:
        for l in lmod:
            for i in range(l):
                inp.write("%d "%f_atom[count+i])
            inp.write("%s\n"%modred_code)
            count += l
        inp.write("\n")

    if 'basis' in lines[10]:
        for l in lines[11:]:
            inp.write("%s"%l)
        inp.write('\n')

    if 'scrf' in lines[8]:
        inp.write('radii=uff\nalpha=1.2\neps=4.0\n\n')
    
    inp.close()


def gen_pdbfiles(wdir,step,tmppdb):
    new_dir = '%s/step%spdbs'%(wdir,step)
    if os.path.isdir(new_dir):
        system_run('rm -r %s'%new_dir)
        system_run('mkdir %s'%new_dir)
    else:
        system_run('mkdir %s'%new_dir)
    os.chdir(new_dir)
    system_run('gopt_to_pdb.py %s %s/step-%s-out 0'%(tmppdb,wdir,step))
    os.chdir('%s'%wdir)
    pdb_name = []
    for pdbf in glob('%s/*.pdb'%new_dir):
        m = re.search('(\d+).pdb',pdbf)
        pdb_name.append( int(m.group(1)) )
    return max(pdb_name), new_dir


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Prepare template PDB files, write input files, save output PDB files in working directory')
    parser.add_argument('-step', dest='step', default=0, type=int, 
    help='step0: read noh, addh pdbs, write_final_pdb and read input_template write_first_inp,\
          step1: read outputwrite_modred_inp, input_template, write_second_inp,\
          step2: read outputwrite_modred_inp, input_template, write_new_inp')
    parser.add_argument('-wdir', dest='output_dir', default=os.path.abspath('./'), help='working dir')
    parser.add_argument('-tmp', dest='tmp_pdb', default=None, help='template_pdb_file')
    parser.add_argument('-noh', dest='no_h_pdb', default=None, help='trimmed_pdb_file')
    parser.add_argument('-adh', dest='h_add_pdb', default=None, help='hadded_pdb_file')
    parser.add_argument('-new', dest='new_pdb', default=None, help='new_pdb_file')
    parser.add_argument('-intmp', dest='input_tmp', default=None, help='template_for_write_input')
    parser.add_argument('-outf', dest='gau_out', default='1.out', help='output_name')
    parser.add_argument('-inpn', dest='inp_name', default='1.inp', help='input_name')
    parser.add_argument('-m', dest='multiplicity', default=1, type=int, help='multiplicity')
    parser.add_argument('-c', dest='lig_charge', default=0, type=int, help='charge_of_ligand')

    args = parser.parse_args()

    step = args.step
    wdir = args.output_dir
    if args.tmp_pdb is None:
        tmp_pdb = '%s/template.pdb'%wdir
    else:
        tmp_pdb   = args.tmp_pdb

#    i_name = []
#    for gau_input in glob('%s/*'%wdir):
#        m = re.search(r'-(\d+)-inp', gau_input)
#        if m:
#            i_name.append( int(m.group(1)) )
#    if len(i_name) > 0:
#        step = max(i_name)
#        if filecmp.cmp('1.inp','%s/step-%s-inp'%(wdir,step)) is False and filecmp.cmp('1.out','%s/step-%s-out'%(wdir,step)) is False:
#            system_run( 'cp 1.inp step-%s-inp'%(step+1) )
#            system_run( 'cp 1.out step-%s-out'%(step+1) )
#            system_run( 'cp 1.chk step-%s-chk'%(step+1) )
#        else:
#            print "check if the files are propagated correctly!"
#        for n in sorted(i_name):
#            pdb_file = gen_pdbfiles(wdir,n,tmppdb)
            
    nohpdb   = args.no_h_pdb
    adhpdb   = args.h_add_pdb
#    newpdb   = args.new_pdb
    int_tmp  = args.input_tmp
    gauout   = args.gau_out
    inp_name = args.inp_name
    multi    = args.multiplicity
    charge   = args.lig_charge
    wdir     = args.output_dir

    if step == 0:
        pic_atom, tot_charge = pdb_after_addh(nohpdb,adhpdb)
        write_pdb('%s'%tmp_pdb,pic_atom)
    elif step == 1:
        i_name = []
        for gau_input in glob('%s/*inp'%wdir):
            m = re.search(r'-(\d+)-inp', gau_input)
            if m:
                i_name.append( int(m.group(1)) )
        if len(i_name) == 0:
            i_step = 1
            system_run( 'cp 1.inp step-%s-inp'%(i_step) )
            system_run( 'cp 1.out step-%s-out'%(i_step) )
            system_run( 'cp 1.chk step-%s-chk'%(i_step) )
        else:
            i_step = max(i_name)
            if filecmp.cmp('1.inp','%s/step-%s-inp'%(wdir,i_step)) is False and filecmp.cmp('1.out','%s/step-%s-out'%(wdir,i_step)) is False:
                i_step += 1
                system_run( 'cp 1.inp step-%s-inp'%(i_step) )
                system_run( 'cp 1.out step-%s-out'%(i_step) )
                system_run( 'cp 1.chk step-%s-chk'%(i_step) )
            elif filecmp.cmp('1.inp','%s/step-%s-inp'%(wdir,i_step)) is True and filecmp.cmp('1.out','%s/step-%s-out'%(wdir,i_step)) is True:
                i_step = i_step
            else:
                print "check if the files are propagated correctly!"
                sys.exit()
#            for n in sorted(i_name):
#                pdb_file, new_dir = gen_pdbfiles(wdir,n,tmp_pdb)
        pdb_file, new_dir = gen_pdbfiles(wdir,i_step,tmp_pdb)

        pic_atom, res_info, tot_charge = read_pdb('%s/%s.pdb'%(new_dir,pdb_file))
    elif step == 2:
        newpdb = args.new_pdb
        pic_atom, res_info, tot_charge = read_pdb(newpdb)
    write_input('%s/%s'%(wdir,inp_name),int_tmp,charge,multi,pic_atom,tot_charge)
        
