#!/usr/bin/env python

from rms import *
from numpy import *
import sys, re, os
import glob

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

def read_gaussian_opt(gfile,natoms,step):
    with open(gfile) as f:
        lines = f.readlines()
    p_start = []
    for i in range(len(lines)):
#        if 'Leave Link 9999' in lines[i]: break
        if 'Standard orientation' in lines[i]:
            p_start.append(i+5)

    ### Read in stardard orientation geom
    opt = {}
    for j in range(len(p_start)):
        key = step+j
        opt[key] = []
        for i in range(p_start[j],p_start[j]+natoms):
            v = lines[i].split()
            opt[key].append([float(v[3]),float(v[4]),float(v[5])])
        opt[key] = array(opt[key])

    return opt


def update_xyz(pdb,xyz):
    sel_atom = []
    for i in range(len(pdb)):
        atom = pdb[i]
        atom[8:11] = xyz[i,:]
        sel_atom.append(atom)
    return sel_atom

def get_ca(pdb):
    map = []
    xyz_i = []
    for i in range(len(pdb)):
        if pdb[i][2].strip() == 'CA':
            map.append(i)
            xyz_i.append(pdb[i][8:11])
    return map, array(xyz_i)

if __name__ == '__main__':
    pdbf = sys.argv[1]
    steps = []
    files = []
    for i in range(2,len(sys.argv),2):
        files.append(sys.argv[i])
        steps.append(int(sys.argv[i+1]))

    pdb, res_info, tot_charge = read_pdb(pdbf)
    map, xyz_i = get_ca(pdb)

    natoms = len(pdb)
    rot_opt = {}
    count = 0
    for k in range(len(files)):
        output = files[k]
        count += int(steps[k])
        step = count
        opt = read_gaussian_opt(output,natoms,step)
        keys = sorted(opt.keys())
        if step in rot_opt.keys() and not array_equal(abs(opt[keys[0]]),abs(rot_opt[step])):
            print "Something is wrong"
            print opt[keys[0]]
            print -1*rot_opt[step]
        for key in sorted(opt.keys()):
            if key not in rot_opt.keys():
                rot_opt[key] = opt[key]
    for key in sorted(rot_opt.keys()):
        xyz_c = rot_opt[key]
        (c_trans,U,ref_trans) = rms_fit(xyz_i,xyz_c[map])
        xyz_n = dot( xyz_c-c_trans, U ) + ref_trans
        sel_atom = update_xyz(pdb,xyz_n)
        name = str(key)+'.pdb'
        write_pdb(name,sel_atom)
