#!/usr/bin/env python

import os, sys, re
from numpy import *
from read_probe import *
from glob import glob


### /home/ndyonker/chem/comt/model-trang/comt-m2-ts-001 ###
### /home/ndyonker/chem/comt/make_pdb ###

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
            charge = line[78:80]
        except:
            charge = '0 '
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
    return pdb, res_info



def gen_pdbfiles(dirf,project_dir,tmppdb):
    new_dir = '%s/pdbfiles'%dirf
    if os.path.isdir(new_dir):
        system_run('rm -r %s'%new_dir)
        system_run('mkdir %s'%new_dir)
    else:
        system_run('mkdir %s'%new_dir)
    os.chdir(new_dir)
    system_run('gopt_to_pdb.py %s %s/1.out 0'%(tmppdb,dirf))
    os.chdir('%s'%project_dir)
    pdb_name = []
    for pdbf in glob('%s/*.pdb'%new_dir):
        m = re.search('(\d+).pdb',pdbf.split('/')[-1])
        pdb_name.append( int(m.group(1)) )
    return max(pdb_name)


if __name__ == '__main__':

    if len(sys.argv) < 3:
        print "Usage write_input.py final_i.pdb project_dir"
        exit()
    tmppdb = sys.argv[1]
    count = int(tmppdb[:-4].split('_')[1])
    project_dir = sys.argv[2]
    dirf = '%s/model%s-ts-001'%(project_dir,count)
    with open('%s/1.inp'%dirf) as f:
        lines = f.readlines()
    charge = lines[7].split()[0]
    
    pdb_name = gen_pdbfiles(dirf,project_dir,tmppdb)
    newpdb = '%s/pdbfiles/%d.pdb'%(dirf,pdb_name)

    tmp_pdb, res_info = read_pdb(tmppdb)
    new_pdb, binfo = read_pdb(newpdb)     #can be just xyz files from cerius or pymol
    #Pic_atom = []
    #For i in range(len(new_pdb)):
    #    line = tmp_xyz[i]
    #    pic_atom.append([line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],new_pdb[i][8],new_pdb[i][9],new_pdb[i][10],line[11],line[12],line[13],' H',line[15],line[16]])

    ### Need to add more type of calculation basis sets etc. ###

    work_dir = '%s/opt'%dirf
    system_run('mkdir %s'%work_dir)
    os.chdir('%s'%work_dir)
    input = open('1.inp','w')
    input.write("%chk=1.chk\n")     #write check file into 1.chk
    input.write("%nproc=20\n")      #use 24 processors
    input.write("%mem=14GB\n")      #use memory
    
    input.write("#P b3lyp/genecp opt(ts,noeigen,calcfc) freq scf=(xqc,maxcon=128,maxcyc=128)\n")      #many things can be changed to user input
    input.write("\n")
    input.write("%d\n"%count)
    input.write("\n")
    input.write("%s 1\n"%charge)        #charge and multiplicity, also need to be changed case by case
    
    
    for i in range(len(new_pdb)):
        line = tmp_pdb[i]
        xyz = new_pdb[i]
        input.write("%4s %6s %8.3f %8.3f %8.3f\n"%(line[14],line[16],xyz[8],xyz[9],xyz[10])) 
    
    input.write("\n")
    input.write("Mg     0                              \n")
    input.write("S   2   1.00                          \n")
    input.write("      0.7250000             -0.4058454\n")
    input.write("      0.1112000              1.1688704\n")
    input.write("S   1   1.00                          \n")
    input.write("      0.0404000              1.0000000\n")
    input.write("P   2   1.00                          \n")
    input.write("      1.2400000             -0.0749753\n")
    input.write("      0.1346000              1.0178183\n")
    input.write("P   1   1.00                          \n")
    input.write("      0.0422000              1.0000000\n")
    input.write("****                                  \n")
    input.write("S     0                               \n")
    input.write("6-31G(d')                             \n")
    input.write("****                                  \n")
    input.write("O     0                               \n")
    input.write("6-31G(d')                             \n")
    input.write("****                                  \n")
    input.write("N     0                               \n")
    input.write("6-31G(d')                             \n")
    input.write("****                                  \n")
    input.write("C     0                               \n")
    input.write("6-31G                                 \n")
    input.write("****                                  \n")
    input.write("H     0                               \n")
    input.write("6-31G                                 \n")
    input.write("****                                  \n")
    input.write("                                      \n")
    input.write("Mg                                    \n")
    input.write("lanl2dz                               \n")

    input.write("\n")
    input.close()
    os.chdir('%s'%project_dir)
