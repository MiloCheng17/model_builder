#!/usr/bin/env python

import os, sys, re
from numpy import *
from glob import glob
import operator


def write_uniquef(dir):
    outfs = sorted(glob('%s/*_rin.dat'%dir))
    for outf in sorted(glob('%s/*_rin.dat'%dir)):
        print outf
        d = genfromtxt(outf,delimiter='>',dtype='str',filling_values='0')
        try:
            dats = vstack((dats,d))
        except:
            dats = d 
    id_mat = dats[:,-1].argsort()
    dats = dats[id_mat]
    
    f1 = open('unique_res_set.dat','w')

    uniq, idx, counts = unique(dats[:,-2:],axis=0,return_index=True,return_counts=True)
    for i in range(len(idx)):
        f1.write('%10d: %10d:'%(idx[i],counts[i]))
        f1.write('%s\n'%dats[idx[i],-1])
    f1.close()

def write_freq_res(uniquef):
    f2 = open('freq_per_res.dat','w')
    res = {}
    with open(uniquef) as f:
        lines = f.readlines()

        for line in lines:
            c = line.split(':')
            for i in map(int,c[-1].split(',')[:-1]):
                if i not in res.keys():
                    res[i] = 1
                else:
                    res[i] += 1
    resf = sorted(res.items(), key=operator.itemgetter(1),reverse=True)
    for i in range(len(resf)):
        f2.write('%6d %8d\n'%(resf[i][0],resf[i][1]))
    f2.close()

        
if __name__ == '__main__':
    dir = sys.argv[1]
    write_uniquef(dir)
    write_freq_res('unique_res_set.dat')

