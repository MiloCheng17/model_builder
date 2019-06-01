#!/usr/bin/env python
"""
This is a program written by Qianyi Cheng in DeYonker Research Group
at University of Memphis.
"""

import os, sys, re


aa_trans_dic = {
'ALA': 'A', 'VAL': 'V', 'PHE': 'F', 'PRO': 'P', 'MET': 'M', 'ILE': 'I', 'LEU': 'L', 'ASP': 'D', 'GLU': 'E', 'LYS': 'K',
'ARG': 'R', 'SER': 'S', 'THR': 'T', 'TYR': 'Y', 'HIS': 'H', 'CYS': 'C', 'ASN': 'N', 'GLN': 'Q', 'TRP': 'W', 'GLY': 'G',
'2AS': 'D', '3AH': 'H', '5HP': 'E', 'ACL': 'R', 'AIB': 'A', 'ALM': 'A', 'ALO': 'T', 'ALY': 'K', 'ARM': 'R', 'ASA': 'D',
'ASB': 'D', 'ASK': 'D', 'ASL': 'D', 'ASQ': 'D', 'AYA': 'A', 'BCS': 'C', 'BHD': 'D', 'BMT': 'T', 'BNN': 'A', 'BUC': 'C',
'BUG': 'L', 'C5C': 'C', 'C6C': 'C', 'CCS': 'C', 'CEA': 'C', 'CHG': 'A', 'CLE': 'L', 'CME': 'C', 'CSD': 'A', 'CSO': 'C',
'CSP': 'C', 'CSS': 'C', 'CSW': 'C', 'CXM': 'M', 'CY1': 'C', 'CY3': 'C', 'CYG': 'C', 'CYM': 'C', 'CYQ': 'C', 'DAH': 'F',
'DAL': 'A', 'DAR': 'R', 'DAS': 'D', 'DCY': 'C', 'DGL': 'E', 'DGN': 'Q', 'DHA': 'A', 'DHI': 'H', 'DIL': 'I', 'DIV': 'V',
'DLE': 'L', 'DLY': 'K', 'DNP': 'A', 'DPN': 'F', 'DPR': 'P', 'DSN': 'S', 'DSP': 'D', 'DTH': 'T', 'DTR': 'W', 'DTY': 'Y',
'DVA': 'V', 'EFC': 'C', 'FLA': 'A', 'FME': 'M', 'GGL': 'E', 'GLZ': 'G', 'GMA': 'E', 'GSC': 'G', 'HAC': 'A', 'HAR': 'R',
'HIC': 'H', 'HIP': 'H', 'HMR': 'R', 'HPQ': 'F', 'HTR': 'W', 'HYP': 'P', 'IIL': 'I', 'IYR': 'Y', 'KCX': 'K', 'LLP': 'K',
'LLY': 'K', 'LTR': 'W', 'LYM': 'K', 'LYZ': 'K', 'MAA': 'A', 'MEN': 'N', 'MHS': 'H', 'MIS': 'S', 'MLE': 'L', 'MPQ': 'G',
'MSA': 'G', 'MSE': 'M', 'MVA': 'V', 'NEM': 'H', 'NEP': 'H', 'NLE': 'L', 'NLN': 'L', 'NLP': 'L', 'NMC': 'G', 'OAS': 'S',
'OCS': 'C', 'OMT': 'M', 'PAQ': 'Y', 'PCA': 'E', 'PEC': 'C', 'PHI': 'F', 'PHL': 'F', 'PR3': 'C', 'PRR': 'A', 'PTR': 'Y',
'SAC': 'S', 'SAR': 'G', 'SCH': 'C', 'SCS': 'C', 'SCY': 'C', 'SEL': 'S', 'SEP': 'S', 'SET': 'S', 'SHC': 'C', 'SHR': 'K',
'SOC': 'C', 'STY': 'Y', 'SVA': 'S', 'TIH': 'A', 'TPL': 'W', 'TPO': 'T', 'TPQ': 'A', 'TRG': 'K', 'TRO': 'W', 'TYB': 'Y',
'TYQ': 'Y', 'TYS': 'Y', 'TYY': 'Y', 'AGM': 'R', 'GL3': 'G', 'SMC': 'C', 'ASX': 'B', 'CGU': 'E', 'CSX': 'C', 'GLX': 'Z',
'O': 'W'}

mc_atoms_dic = {'N': '', 'CA': '', 'C': '', 'O': '', 'H': '', 'HA': '', 'OXT': '', 'HA2': '', 'HA3': '', 'H?': '', 'W': ''}

def get_res_type(resID, atom):
    if resID in aa_trans_dic.keys():
        if atom in mc_atoms_dic.keys():
            side = 'mc'
        else:
            side = 'sc'
    else:
        if resID == 'HOH':
            side = 'solvent'
        else:
            side = 'ligand'
    return side

def get_inttype(c):
    # 'wc':wide contact,'cc': close contact,'so':small overlap
    # 'bo':big overlap
    # 'hb':hydrogen bond
    if c in ['wc','cc']:
        action = 'cnt'
    elif c in ['bo','so']:
        action = 'ovl'
    elif c == 'hb':
        action = 'hbond'
    else:
        print 'Cannot find interaction type!'
        return "None"
    return action

def set_mc_sc_ligand(side1, side2):
    if side1==None or side2==None:
        return None, None
    elif side1 == 'mc':
        return side1, side2
    elif side1 == 'sc':
         if side2 == 'mc':
             return side2, side1
         else:
             return side1, side2
    elif side1 == 'ligand':
        if (side2 in ['mc', 'sc']):
            return side2, side1
        else:
            return side1, side2
    elif side1 == 'solvent':
        return side2, side1

def get_side(c):
    cha    = c[:2].strip()          # chain 
    res_id = c[2:6].strip()              # res_id 
    res_nm = c[6:10].strip()     # res_name
    atom   = c[10:-1].strip()      # atom_name
    if c[-1] == ' ':
        spe_c = 'A'
    else:
        spe_c = c[-1]
    return cha, res_id, res_nm, atom, spe_c

def probe2Sif(probefile):
    res_parts = {}
    res_parts1 = {}

    with open(probefile,'r') as f:
        lines = f.readlines()

    for line in lines:
        c = line.split(':')
        obj1 = c[1]
        acts = c[2]

        cha1, res_id1, res_nm1, atom1, spe_c1 = get_side(c[3])           
        cha2, res_id2, res_nm2, atom2, spe_c2 = get_side(c[4])           

        if cha1 == cha2 and res_id1 == res_id2: continue
        ### for generate sif file ###
        part1 = cha1+':'+res_id1+':'+'_'+':'+res_nm1
        part2 = cha2+':'+res_id2+':'+'_'+':'+res_nm2
        side1 = get_res_type(res_nm1,atom1)
        side2 = get_res_type(res_nm2,atom2)
        action = get_inttype(acts)
        sideL, sideR = set_mc_sc_ligand(side1,side2)
        act = action+':'+sideL+'_'+sideR

        if  (part1,act,part2) not in res_parts.keys() and (part2,act,part1) in res_parts.keys():
            key = (part2,act,part1)
            if side1 not in res_parts[key][1]:
                res_parts[key][1].append(side1)
            if side2 not in res_parts[key][0]:
                res_parts[key][0].append(side2)
            key1 = (part2,part1)

        elif (part1,act,part2) not in res_parts.keys() and (part2,act,part1) not in res_parts.keys():
            key = (part1,act,part2)
            res_parts[key] = [[side1],[side2]]
            key1 = (part1,part2)

        else:
            if side1 not in res_parts[key][0]:
                res_parts[key][0].append(side1)
            if side2 not in res_parts[key][1]:
                res_parts[key][1].append(side2)
            key1 = (part1,part2)

    interaction = {}
    combi_interaction = []

    for key in res_parts.keys():
        if key not in interaction:
             interaction[key] = res_parts[key]

        k1 = [key[0], key[2]]
        k2 = [key[2], key[0]]
        if (k1 not in combi_interaction) and (k2 not in combi_interaction):
            key1 = (key[0],'combi:all_all',key[2])
            interaction[key1] = []
            combi_interaction.append(k1)
    
    f = open(probefile.replace(".probe", ".sif"), "w")
    for key in interaction.keys():
        if ('combi' not in key[1]):
            act = key[0]+' '+key[1]+' '+key[2]
            f.write(act+"\n")

    for key in interaction.keys():
        if ('combi' in key[1]):
            act = key[0]+' '+key[1]+' '+key[2]
            f.write(act+"\n")

    f.close()        

if __name__ == "__main__":
    if (len(sys.argv) < 2):
        print "Usage: probe_file"
        exit()

    probe2Sif(sys.argv[1])

