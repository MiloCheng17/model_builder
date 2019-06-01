#!/usr/bin/env python
"""
This is a program written by Qianyi Cheng in DeYonker Research Group
at University of Memphis.
"""

from rms import *
from numpy import *
import sys, re, os
from read_write_pdb import *
import argparse

def system_run(cmd):
    print cmd
    exit = os.system(cmd)
    if ( exit != 0 ):
        print 'failed to run:'
        print cmd
        sys.exit()
 
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

if __name__ == '__main__':
    """ Usage: gopt_to_pdb.py -o ../1.out -s -1 -p ../template.pdb """
    parser = argparse.ArgumentParser(description='generate pdbfiles from 1.out')
    parser.add_argument('-f', dest='frame',default=None,help='select frame/range')
    parser.add_argument('-o', dest='outputs',default='1.out',help='output files')
    parser.add_argument('-s', dest='start',default=0,type=int,help='output steps')
    parser.add_argument('-p', dest='pdbf',default='template.pdb',help='template pdb file')

    args = parser.parse_args()
    steps = []
    files = []
    if ',' in args.outputs:
        files = args.outputs.split(',')
        steps = args.start.split(',')
    else:
        files = [args.outputs]
        steps = [args.start]
    pdbf = args.pdbf
#    pdbf = sys.argv[1]
#    steps = []
#    files = []
#    for i in range(2,len(sys.argv),2):
#        files.append(sys.argv[i])
#        steps.append(int(sys.argv[i+1]))

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
    if args.frame is None:
        for key in sorted(rot_opt.keys()):
            xyz_c = rot_opt[key]
            (c_trans,U,ref_trans) = rms_fit(xyz_i,xyz_c[map])
            xyz_n = dot( xyz_c-c_trans, U ) + ref_trans
            sel_atom = update_xyz(pdb,xyz_n)
            name = str(key)+'.pdb'
            write_pdb(name,sel_atom)
    elif args.frame == '-1':
        key = sorted(rot_opt.keys())[-1]
        xyz_c = rot_opt[key]
        (c_trans,U,ref_trans) = rms_fit(xyz_i,xyz_c[map])
        xyz_n = dot( xyz_c-c_trans, U ) + ref_trans
        sel_atom = update_xyz(pdb,xyz_n)
        name = str(key)+'.pdb'
        write_pdb(name,sel_atom)
    else:
        print "The frame is not clear!"
        sys.exit()
