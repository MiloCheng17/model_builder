#!/usr/bin/env python

from read_write_pdb import *
from rms import *
from numpy import *
import sys, re, os
import argparse

def get_gaussian_freq(gfile):
    with open(gfile) as f:
        lines = f.readlines()
    p_start = []
    f_lines = []
    atom_name = []
    for i in range(len(lines)):
        if 'NAtoms=' in lines[i]:
            num_atoms = int(lines[i].split()[1])
            break
    for i in range(len(lines)):
        if 'Charge =' in lines[i]:
            v = lines[i].split()
            charge = int(v[2])
            multip = int(v[5])
            if "Redundant internal coordinates found in file." in lines[i+1]:
                for idx in range(i+2,i+2+num_atoms):
                    if ',' in lines[idx]:
                        atom_name.append(lines[idx][:-1].split(','))
                    else:
                        va = lines[idx][:-1].split()
                        if len(va) == 4:
                            atom_name.append([va[0],'0'])
                        else:
                            atom_name.append(lines[idx][:-1].split()[:2])

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
    print len(opt)

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

    return opt, atom_idx[3], freq_xyz, freq_info, charge, multip, atom_name        


def form_irc_xyz(opt,atom_idx,xyz,scale):
    irc1 = []
    irc2 = []
#    scale = 0.1
    for i in range(opt.shape[0]):
        if i not in atom_idx:
            irc1.append(opt[i,:])
            irc2.append(opt[i,:])
        else:
            irc1.append(opt[i,:]+scale*array(xyz[atom_idx.index(i)]))
            irc2.append(opt[i,:]-scale*array(xyz[atom_idx.index(i)]))
    return array(irc1), array(irc2)

def write_input(inp_f,dir1,dir2,irc1,irc2,atom_name,charge,multip):
    with open(inp_f) as f:
        lines = f.readlines()
    f1 = open('%s/1.inp'%dir1,'w')
    f2 = open('%s/1.inp'%dir2,'w')
    
    for i in range(len(lines)):
        v = lines[i].split()
#        if len(v) == 5 and '#P' not in lines[i] and i > 7:
        if len(v) != 0 and "#P" in v[0]:
#            print v
            for tag in v:
                if 'opt' in tag: 
                    v[v.index(tag)] = 'opt'
                if 'geom' in tag or 'guess' in tag or 'iop' in tag:
                    v[v.index(tag)] = ''
            for word in v:
                f1.write('%s '%word)
                f2.write('%s '%word)
            f1.write('\n')
            f2.write('\n')
        elif len(v) == 2 and v[0] == str(charge) and v[1] == str(multip):
            f1.write(lines[i])
            f2.write(lines[i])
            for atom in range(len(irc1)):
                f1.write('%6s%6s%12.6f%12.6f%12.6f\n'%(atom_name[atom][0],atom_name[atom][1],irc1[atom,0],irc1[atom,1],irc1[atom,2]))
                f2.write('%6s%6s%12.6f%12.6f%12.6f\n'%(atom_name[atom][0],atom_name[atom][1],irc2[atom,0],irc2[atom,1],irc2[atom,2]))
        elif len(v) == 4 or len(v) == 5:
            if len(v[-1]) != 1:
                continue
            else:
                f1.write(lines[i])
                f2.write(lines[i])
        else:
            f1.write(lines[i])
            f2.write(lines[i])
    f1.close()
    f2.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Write irc inputs')
    parser.add_argument('-wdir', dest='output_dir', default=os.path.abspath('./'), help='working dir')
    parser.add_argument('-dir1', dest='dir1', default=os.path.abspath('./irc1'), help='irc1 dir')
    parser.add_argument('-dir2', dest='dir2', default=os.path.abspath('./irc2'), help='irc2 dir')
    parser.add_argument('-outf', dest='gau_out', default=None, help='output_name')
    parser.add_argument('-inpf', dest='gau_inp', default=None, help='input_name')
    parser.add_argument('-s', dest='scale', type=float,default=0.1, help='scale_factor')
    parser.add_argument('-n', dest='num_freq', type=int,default=1, help='scale_factor')

    args = parser.parse_args()
    wdir = args.output_dir
    scale = args.scale
    num_freq = args.num_freq
    dir1 = args.dir1
    dir2 = args.dir2
    if args.gau_out is None:
        out_f = '%s/1.out'%wdir
    else:
        out_f = args.gau_out

    if args.gau_inp is None:
        inp_f = '%s/1.inp'%wdir
    else:
        inp_f = args.gau_inp

    opt, atom_idx, freq_xyz, freq_info, charge, multip, atom_name = get_gaussian_freq(out_f)
#    print charge, multip

    xyz = freq_xyz[num_freq]
    irc1, irc2 = form_irc_xyz(opt,atom_idx,xyz,scale)
    write_input(inp_f,'irc1','irc2',irc1,irc2,atom_name,charge,multip)
