#!/usr/bin/env python

from PDB import *
from rms import *
from numpy import *
import sys, re, os
import glob

def get_gaussian_freq(gfile,num_atoms):
    with open(gfile) as f:
        lines = f.readlines()
    p_start = []
    f_lines = []
    for i in range(len(lines)):
        if 'Standard orientation' in lines[i]:
            p_start.append(i+5)
        if 'Frequencies' in lines[i][1:12]:
            f_lines.append(i)

    ### Read in stardard orientation
    opt = []
    for i in range(p_start[-1],p_start[-1]+num_atoms):
        v = lines[i].split()
        opt.append([float(v[3]),float(v[4]),float(v[5])])
    opt = array(opt)

    ### Read in imaginary frequency
    freq_xyz  = {}
    freq_info = {}
    keys = []
    num_atoms = f_lines[1]-f_lines[0]-7
    atom_idx = {}
    for i in range(len(f_lines)):
        start = f_lines[i] - 2
        end = f_lines[i]+4+num_atoms
        modes = lines[start].split()
        sym   = lines[start+1].split()
        mag   = lines[start+2].split()[2:]
        red   = lines[start+3].split()[3:]
        force = lines[start+4].split()[3:]
        ir_it = lines[start+5].split()[3:]
        for j in range(3):
            key = int(modes[j])
            keys.append(key)
            freq_info[key] = (sym[j],float(mag[j]),float(red[j]),float(force[j]),float(ir_it[j])) 
            freq_xyz[key]  = []
        atom_idx[keys[-1]] = []
        for j in range(start+7,end+1):
            v = lines[j].split()
            atom_idx[keys[-1]].append(int(v[0])-1)
            for k in range(3):
                freq_xyz[keys[i*3+k]].append([float(v[3*k+2]),float(v[3*k+3]),float(v[3*k+4])])

    return opt, atom_idx[3], freq_xyz, freq_info        


def form_irc_xyz(opt,atom_idx,xyz):
    irc1 = []
    irc2 = []
#    scale = 0.5
    scale = 0.1
    for i in range(opt.shape[0]):
        if i not in atom_idx:
            irc1.append(opt[i,:])
            irc2.append(opt[i,:])
        else:
            irc1.append(opt[i,:]+scale*array(xyz[atom_idx.index(i)]))
            irc2.append(opt[i,:]-scale*array(xyz[atom_idx.index(i)]))
    return array(irc1), array(irc2)

def write_input(input,dir1,dir2,irc1,irc2):
    with open(input) as f:
        lines = f.readlines()
    f1 = open('%s/1.inp'%dir1,'w')
    f2 = open('%s/1.inp'%dir2,'w')
    count = 0
    for i in range(len(lines)):
        v = lines[i].split()
        if len(v) == 5 and '#P' not in lines[i] and i > 7:
            f1.write('%6s%6s%12.6f%12.6f%12.6f\n'%(v[0],v[1],irc1[count,0],irc1[count,1],irc1[count,2]))
            f2.write('%6s%6s%12.6f%12.6f%12.6f\n'%(v[0],v[1],irc2[count,0],irc2[count,1],irc2[count,2]))
            count += 1
        else: 
            f1.write(lines[i])
            f2.write(lines[i])
    f1.close()
    f2.close()

def update_xyz(pdb,xyz):
    sel_atom = []
    for i in range(len(pdb)):
        atom = pdb[i]
        atom[8:11] = xyz[i,:]
        sel_atom.append(atom)
    return sel_atom


if __name__ == '__main__':
    outdir = sys.argv[1]
    num_atoms=int(sys.argv[2])
    input = sys.argv[3]
    opt, atom_idx, freq_xyz, freq_info = get_gaussian_freq(outdir,num_atoms)

    xyz = freq_xyz[1]
    irc1, irc2 = form_irc_xyz(opt,atom_idx,xyz)
    write_input(input,'irc1','irc2',irc1,irc2)
