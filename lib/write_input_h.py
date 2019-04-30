#!/usr/bin/env python

import os, sys, re
from numpy import *
from read_probe import *


### /home/ndyonker/chem/comt/model-trang/comt-m2-ts-001 ###
### /home/ndyonker/chem/comt/make_pdb ###


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
        #print charge, charge.isdigit()
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
    return pdb, res_info, tot_charge    ### suppose pymol charge is correct



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


def write_input(tmp_input,charge,pic_atom,tot_charge):
    with open(tmp_input) as f:
        lines = f.readlines()
    input = open('%d.input'%count,'w')
    input.write("%chk=1.chk\n")     #write check file into 1.chk
    input.write("%nproc=24\n")      #use 24 processors
    input.write("%mem=14GB\n")      #use memory
    
    input.write("#P b3lyp/gen opt freq scf=(xqc,maxcon=128,maxcyc=128)\n")      #many things can be changed to user input
    input.write("\n")
    title = raw_input("Input the descriptive name of this calculation: ")
    input.write("%s\n"%title)
    input.write("\n")
    input.write("%d 1\n"%(charge+tot_charge))        #charge and multiplicity, also need to be changed case by case
    
    
    for i in range(len(pic_atom)):
        input.write("%4s %6s %8.3f %8.3f %8.3f\n"%(atom[i],hold[i],xyz[i][0],xyz[i][1],xyz[i][2])) 
    
    input.write("\n")
    input.write("O N S\n")
    input.write("6-31G(d')\n")
    input.write("****\n")
    input.write("C H\n")
    input.write("6-31G\n")
    input.write("****\n")
    input.write("\n")
    input.close()


if __name__ == '__main__':

    if len(sys.argv) < 3:
        print "Usage ~/projects/model_build/lib-2018/write_input.py template.pdb H-added.pdb"
        exit()
    tmppdb = sys.argv[1]
    count = int(tmppdb[:-4].split('_')[1])
    newpdb = sys.argv[2]
    outf = 'final_%d.pdb'%count

    tmp_pdb, res_info, tot_charge_t = read_pdb(tmppdb)
    tmp_xyz = []
    for i in tmp_pdb:
        tmp_xyz.append([i[8],i[9],i[10]])
    new_pdb, binfo, tot_charge = read_pdb(newpdb)     #can be just xyz files from cerius or pymol
    print "total charge is %d" %tot_charge
    pic_atom = []
    xyz = []
    atom = []
    hold = []
    for line in new_pdb:
        xyz.append([line[8],line[9],line[10]])
        if [line[8],line[9],line[10]] not in tmp_xyz:
            atom.append(' H')
            pic_atom.append([line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],' H',line[15],line[16]])
            #pic_atom.append([line[0],line[1],'H',line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],' H',line[15],line[16]])
            hold.append(0)
        else:
            idx = tmp_xyz.index([line[8],line[9],line[10]])
            #line[14] = tmp_pdb[idx][14]
            line = tmp_pdb[idx]
            atom.append(line[14])
            if [line[2],line[5],line[6]] in res_info:
                pic_atom.append([line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],'-1'])
                hold.append(-1)
            else:
                pic_atom.append([line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],line[16]])
                hold.append(0)
    
    write_pdb(outf,pic_atom)

    ### Need to add more type of calculation basis sets etc. ###

    input = open('%d.input'%count,'w')
    input.write("%chk=1.chk\n")     #write check file into 1.chk
    input.write("%nproc=10\n")      #use 24 processors
    input.write("%mem=3800MB\n")      #use memory
    
   # input.write("#P b3lyp/genecp opt freq scf=(xqc,maxcon=128,maxcyc=128)\n")      #many things can be changed to user input
    input.write("#P pm6 opt(calcfc) scf=(xqc,maxcon=128,maxcyc=128)\n")      #many things can be changed to user input
#    input.write("#P pm6 opt scf=(xqc,maxcon=128,maxcyc=128)\n")      #many things can be changed to user input
    input.write("\n")
    input.write("%d\n"%count)
    input.write("\n")
    input.write("%d 1\n"%(1+tot_charge))        #charge and multiplicity, also need to be changed case by case
    
    
    for i in range(len(pic_atom)):
        if atom[i] == " H":
            hold[i] = 0
        else:
            hold[i] = -1
        input.write("%4s %6s %8.3f %8.3f %8.3f\n"%(atom[i],hold[i],xyz[i][0],xyz[i][1],xyz[i][2])) 
    
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

#    input.write("O N S\n")
#    input.write("6-31G(d')\n")
#    input.write("****\n")
#    input.write("C H\n")
#    input.write("6-31G\n")
#    input.write("****\n")
    input.write("\n")
    input.close()
