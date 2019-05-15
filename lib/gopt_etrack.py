#!/usr/bin/env python

from numpy import *
import sys, re, os
import glob

def read_gaussian_opt(gfile):
    ### step: 1--n
    ### SCF Done: 1--n+1 (freq also have energy)
    ### Maximum Force: 1--n+1 (freq also have)
    ### Eigenvalues not sure
    ### Standard orientation: 1--n+2 (initial and freq also have)
    with open(gfile) as f:
        lines = f.readlines()
    p_start = []
    steps   = []
    force   = {}
    scf     = {}
    eigen   = {}
    count = 0
    for i in range(len(lines)):
        if 'NAtoms=' in lines[i]:
            natoms = int(lines[i].split()[1])
        if 'Leave Link 9999' in lines[i]: break
        if 'Standard orientation' in lines[i]:
            p_start.append(i+5)
        if 'SCF Done' in lines[i]:
            count += 1
            scf[count] = (float(lines[i].split()[4]))
        if 'out of a maximum of ' in lines[i]:
            stp = int(lines[i].split()[2])
            steps.append(stp)
        if 'ITU= ' in lines[i]:
            eigen[stp] = []
            for j in range(1,5):
                if 'Eigenvalues ---' in lines[i+j]:
                    for eigen_value in lines[i+j].split()[2:]:
                        eigen[stp].append(float(eigen_value))
                    break
        if 'Maximum Force' in lines[i]:
            force[count] = []
            for k in range(4):
                v = lines[i+k].split()
                force[count].append(float(v[2]))
                force[count].append(float(v[3]))
                force[count].append(v[4])

    ### Read in stardard orientation geom
    opt = {}
    for j in range(len(p_start)):
        key = j
        opt[key] = []
        for i in range(p_start[j],p_start[j]+natoms):
            v = lines[i].split()
            opt[key].append([float(v[3]),float(v[4]),float(v[5])])
        opt[key] = array(opt[key])

    return opt, force, scf, eigen

def print_etrack(force, scf, eigen):
    for num in sorted(scf.keys()):
        print 'Step %4d: SCF Done:  E = %16.8f'%(num,scf[num])
        if num in force.keys():
            print '        Item               Value     Threshold  Converged?'
            print 'Maximum Force            %9.6f %11.6f     %-4s'%(force[num][0],force[num][1],force[num][2])
            print 'RMS     Force            %9.6f %11.6f     %-4s'%(force[num][3],force[num][4],force[num][5])
            print 'Maximum Displacement     %9.6f %11.6f     %-4s'%(force[num][6],force[num][7],force[num][8])
            print 'RMS     Displacement     %9.6f %11.6f     %-4s'%(force[num][9],force[num][10],force[num][11])
        else:
            print 'Step %d has no force value yet!'%num
        if num in eigen.keys() and len(eigen[num]) != 0:
            print 'Eigenvalues of step %4d: '%(num), tuple(eigen[num])
        else:
            print 'Step %d has no Eigenvalues yet!'%num



def write_etrack(force, scf, eigen):
    efile = open('etrack.dat','w')
    for num in sorted(scf.keys()):
        efile.write('Step %4d: '%(num))
        efile.write('SCF Done:  E = %16.8f\n'%scf[num])
        if num in force.keys() and len(force[num]) != 0:
            efile.write('        Item               Value     Threshold  Converged?\n')
            efile.write('Maximum Force            %9.6f %11.6f     %-4s\n'%(force[num][0],force[num][1],force[num][2]))
            efile.write('RMS     Force            %9.6f %11.6f     %-4s\n'%(force[num][3],force[num][4],force[num][5]))
            efile.write('Maximum Displacement     %9.6f %11.6f     %-4s\n'%(force[num][6],force[num][7],force[num][8]))
            efile.write('RMS     Displacement     %9.6f %11.6f     %-4s\n'%(force[num][9],force[num][10],force[num][11]))
        if num in eigen.keys() and len(eigen[num]) != 0:
            efile.write('Eigenvalues of step %4d:'%(num))
            for eigen_vl in eigen[num]:
                efile.write('%11.5f'%eigen_vl)
            efile.write('\n')
    efile.close()


if __name__ == '__main__':
    output = sys.argv[1]
    opt, force, scf, eigen = read_gaussian_opt(output)
    print_etrack(force,scf,eigen)
    write_etrack(force,scf,eigen)
