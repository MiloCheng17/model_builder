#!/usr/bin/env python

from PDB import *
from rms import *
from numpy import *
import sys, re, os
import glob

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

    pdb = read_pdb(pdbf)
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
