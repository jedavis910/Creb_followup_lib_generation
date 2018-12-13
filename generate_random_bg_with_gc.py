#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 18:35:18 2018

@author: Leelu
"""
import random
import re

random.seed(123)


def csv_writer(lib, filename):
	with open(filename, 'w') as outfile:
		for x in lib:
			outfile.write(x+','+lib[x]+'\n')

#Design of the background GC content library
    
#Generate random bgs, calculate their GC content
    
def random_bgs(n):
	bgs = []
	for i in range(n):
		bg = ''.join([random.choice(['A', 'C', 'G', 'T']) for j in range(150)])
		bgs.append(bg)
	return bgs

bgs = random_bgs(100000)


def check_REs_bg(lib):
    '''
    Remove any sequences that contain restriction sites used for cloning
    '''
    good_seqs = []
    
    # just check by hand, easier
    # MluI, ACGCGT
    # KpnI, GGTACC
    # XbaI, TCTAGA
    # SpeI, ACTAGT
    # CRE half site, CGTCA and TGACG
    
    patterns = [re.compile('ACGCGT'), re.compile('GGTACC'), 
                re.compile('TCTAGA'), re.compile('ACTAGT'),
                re.compile('CGTCA'), re.compile('TGACG')]
    
    for x in lib:
        # trim sequences so as not to include the KpnI and MluI site that were added in
        matches = [re.search(pattern, x) for pattern in patterns]
        match_count = sum([1 for match in matches if match != None])
        if match_count == 0:
            good_seqs.append(x)
            
    print(len(lib) - len(good_seqs), "sequences removed.")
    print(len(good_seqs), "sequences after removing those with restriction sites.")
        
    return good_seqs

bgs_noRE = check_REs_bg(bgs)


def calculate_gc(seq):
	count = 0
	for i in range(len(seq)):
		if seq[i] == 'G' or seq[i] == 'C':
			count += 1

	return count / float(len(seq))

bgs_noRE_gc = {i : str(calculate_gc(i)) for i in bgs_noRE}


csv_writer(bgs_noRE_gc, 'bgs_noRE_gc.csv')


