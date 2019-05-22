#!/usr/bin/env python

import os, sys, re
from numpy import *

#outf = sys.argv[1]

#scfe = []
def read_gaussian_opt(gfile):
    with open(gfile) as f:
        fline = f.read().replace('\n ','')
    m = re.search('NImag=(\w+)',fline)
    nimag = int(m.group(0).split('=')[-1])
    m = re.search('NAtoms=(\s+)(\d+)',fline)
    natoms = int(m.group(0).split('=')[-1])
    m = re.search('NBasis=(\s+)(\d+)',fline)
    nbasis = int(m.group(0).split('=')[-1])
#    m = re.search("Zero-point correction=(\s+)(\d+).(\d+)",fline)
#    l = m.groups()[1:]
#    if len(l) == 2: 
#        zero = float(l[0]+'.'+l[1])
#    elif len(l) == 3:
#        zero = float(l[0]+l[1]+'.'+l[2]) 
#
#    print zero
##        if "Thermal correction to Energy=" in l:
#    m = re.search("Thermal correction to Energy=(\s+)(\d+).(\d+)",fline)
#    l = m.groups()[1:]
#    if len(l) == 2: 
#        ther = float(l[0]+'.'+l[1])
#    elif len(l) == 3:
#        ther = float(l[0]+l[1]+'.'+l[2]) 
#    print ther
##        if "Thermal correction to Enthalpy=" in l:
#    m = re.search("Thermal correction to Enthalpy=(\s+)(\d+).(\d+)",fline)
#    l = m.groups()[1:]
#    if len(l) == 2: 
#        he = float(l[0]+'.'+l[1])
#    elif len(l) == 3:
#        he = float(l[0]+l[1]+'.'+l[2]) 
#    print he
##        if "Thermal correction to Gibbs Free Energy=" in l:
#    m = re.search("Thermal correction to Gibbs Free Energy=(\s+)(\d+).(\d+)",fline)
#    l = m.groups()[1:]
#    if len(l) == 2: 
#        ge = float(l[0]+'.'+l[1])
#    elif len(l) == 3:
#        ge = float(l[0]+l[1]+'.'+l[2]) 
#    print ge
##        if "E (Thermal)             CV                S" in l:
##    m = re.search("E (Thermal)             CV                S(\s+)?(\d+).(\d+)",fline)
##    l = m.groups()[1:]
##    if len(l) == 2: 
##        se = float(l[0]+'.'+l[1])
##    elif len(l) == 3:
##        se = float(l[0]+l[1]+'.'+l[2]) 
##        if "Sum of electronic and zero-point Energies" in l:
#    m = re.search(r"Sum of electronic and zero-point Energies=(\s+)(\S)(\d+).(\d+)",fline)
#    l = m.groups()[1:]
#    if len(l) == 2: 
#        szero = float(l[0]+'.'+l[1])
#    elif len(l) == 3:
#        szero = float(l[0]+l[1]+'.'+l[2]) 
#    print szero
##        if "Sum of electronic and thermal Energies" in l:
#    m = re.search("Sum of electronic and thermal Energies=(\s+)(\S)(\d+).(\d+)",fline)
#    l = m.groups()[1:]
#    if len(l) == 2: 
#        sele = float(l[0]+'.'+l[1])
#    elif len(l) == 3:
#        sele = float(l[0]+l[1]+'.'+l[2]) 
##        if "Sum of electronic and thermal Enthalpies" in l:
#    m = re.search("Sum of electronic and thermal Enthalpies=(\s+)(\S)(\d+).(\d+)",fline)
#    l = m.groups()[1:]
#    if len(l) == 2: 
#        enthalpy = float(l[0]+'.'+l[1])
#    elif len(l) == 3:
#        enthalpy = float(l[0]+l[1]+'.'+l[2]) 
##        if "Sum of electronic and thermal Free Energies" in l:
#    m = re.search("Sum of electronic and thermal Free Energies=(\s+)(\S)(\d+).(\d+)",fline)
#    l = m.groups()[1:]
#    if len(l) == 2: 
#        freeg = float(l[0]+'.'+l[1])
#    elif len(l) == 3:
#        freeg = float(l[0]+l[1]+'.'+l[2]) 
        

    scfe = []
    with open(gfile) as f:
        lines = f.readlines()
    for i in range(len(lines)):
        l = lines[i]
        if "SCF Done" in l:
            scfe.append(float(l.split()[4]))
        #Zero-point correction=                           2.285998 (Hartree/Particle)
        #Thermal correction to Energy=                    2.407590
        #Thermal correction to Enthalpy=                  2.408534
        #Thermal correction to Gibbs Free Energy=         2.131712
        #Sum of electronic and zero-point Energies=          -6682.107609
        #Sum of electronic and thermal Energies=             -6681.986017
        #Sum of electronic and thermal Enthalpies=           -6681.985073
        #Sum of electronic and thermal Free Energies=        -6682.261895
        #
        #                    E (Thermal)             CV                S
        #                     KCal/Mol        Cal/Mol-Kelvin    Cal/Mol-Kelvin
        #Total                 1510.786            480.805            582.621
        if "Zero-point correction="  in l:
            zero = float(l.split()[-2])
        if "Thermal correction to Energy=" in l:
            ther = float(l.split()[-1])
        if "Thermal correction to Enthalpy=" in l:
            he = float(l.split()[-1])
        if "Thermal correction to Gibbs Free Energy=" in l:
            ge = float(l.split()[-1])
        if "E (Thermal)             CV                S" in l:
            se = float(lines[i+2].split()[-1])
        if "Sum of electronic and zero-point Energies" in l:
            szero = float(l.split()[-1])
        if "Sum of electronic and thermal Energies" in l:
            sele = float(l.split()[-1])
        if "Sum of electronic and thermal Enthalpies" in l:
            enthalpy = float(l.split()[-1])
        if "Sum of electronic and thermal Free Energies" in l:
            freeg = float(l.split()[-1])
    return scfe[-1], zero, ther, he, ge, se, szero, sele, enthalpy, freeg, nimag, nbasis, natoms   

if __name__ ==  '__main__':
#########################################################################
###     Usage: Program output   ####
### or  Usage: Program output 1 #### (print out numbers for table)
### or  Usage: Program output 2 #### (print out my numbers)
#########################################################################

    outf = sys.argv[1]
    scfe, zero, ther, he, ge, se, szero, sele, enthalpy, freeg, nimag, nbasis, natoms = read_gaussian_opt(outf)
    if len(sys.argv) == 2:
        print "%-50s %30.6f"%("SCFDone:E(RB3LYP)= ",scfe)
        print "%-50s %30.6f"%(" Sum of electronic and zero-point Energies= ",szero)
        print "%-50s %30.6f"%(" Sum of electronic and thermal Energies= ",sele)
        print "%-50s %30.6f"%(" Sum of electronic and thermal Enthalpies= ",enthalpy)
        print "%-50s %30.6f"%(" Sum of electronic and thermal Free Energies= ",freeg)
        print "%-50s %30d"%("NImag= ",nimag)
        print "%-50s %30d"%("NBasis= ",nbasis)
        print "%-50s %30d"%("Natoms= ",natoms)
    elif len(sys.argv) == 3 and sys.argv[2] == "1":   
        print "E(RB3LYP)", scfe, szero, sele, enthalpy, freeg, natoms, nbasis, "NImag", nimag
    elif len(sys.argv) == 3 and sys.argv[2] == "2":
        print scfe, zero, ther, he, ge, se, natoms

