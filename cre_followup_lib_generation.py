import re
import random
import numpy as np


random.seed(123)


def fasta_reader(filename):
	'''
	Read in fasta file as dictionary, with headers as keys and 
	sequences as values
	'''
	infile = open(filename)

	seqs = {}

	while 1:
		line = infile.readline()
		if not line:
			break
		# skip first character '>'
		header = line.strip()[1:]
		seq = infile.readline().strip()
		seqs[header] = seq

	return seqs

def csv_reader(filename):
	infile = open(filename)
	seqs = {}

	for line in infile:
		name, seq = line.strip().split(',')
		seqs[name] = seq

	return seqs

def fasta_writer(outfile, seqs):
	for x in seqs:
		outfile.write('>'+x+'\n'+seqs[x]+'\n')
	outfile.close()

def check_REs(lib):
	'''
	Remove any sequences that contain restriction sites used for cloning
	'''
	good_seqs = {}

	# just check by hand, easier
	# MluI, ACGCGT
	# KpnI, GGTACC
	# XbaI, TCTAGA
	# SpeI, ACTAGT

	patterns = [re.compile('ACGCGT'), re.compile('GGTACC'), 
		re.compile('TCTAGA'), re.compile('ACTAGT')]


	for x in lib:
		# trim sequences so as not to include the KpnI and MluI site that were added in
		matches = [re.search(pattern, lib[x][25:-25]) for pattern in patterns]
		match_count = sum([1 for match in matches if match != None])
		if match_count == 0:
			good_seqs[x] = lib[x]
		else:
			print(x, "removed")

	print(len(lib) - len(good_seqs), "sequences removed.")
	print(len(good_seqs), "sequences after removing those with restriction sites.")
	
	return good_seqs

def reverse_complement(seq):
	"""
	Return the reverse complement of a nucleotide string
	"""
	complement = {'A': 'T', 'T':'A', 'C':'G', 'G':'C'}
	rc = ''.join([complement[nt] for nt in seq[::-1]])
	
	return rc


def best_A_content(oligo):
	'''
	Choose the strand with the lowest A content because A's are harder to
	synthesize.
	'''

	# get the reverse compliment of the oligo
	rc_oligo = reverse_complement(oligo)

	# count the number of A's for both strands
	oligo_As = sum( [1 for nt in oligo if nt == 'A'] )
	rc_As = sum( [1 for nt in rc_oligo if nt == 'A'] )

	# save the oligo with lower number of A's
	if oligo_As < rc_As:
		final_oligo = oligo
	else:
		final_oligo = rc_oligo

	return final_oligo

def csv_writer(lib, filename):
	with open(filename, 'w') as outfile:
		for x in lib:
			outfile.write(x+','+lib[x]+'\n')
            
def capitalize_values(d):
    result = {}
    for key, value in d.items():
        upper_value = value.upper()
        result[key] = upper_value
    return result

def similar(a, b):
    # must be same length sequence
    num_same = sum([1 for i in range(len(a)) if a[i] == b[i]])
    return num_same

def most_scramble(sequence, n_iter):
	'''
	Scramble sequence n_iter times and choose the most scrambled sequence
	'''
	scrambled_seqs = [''.join(random.sample(list(sequence), 
                                         len(sequence))) for i in range(n_iter)]
	similarity = [similar(sequence, x) for x in scrambled_seqs]
	max_scrambled = scrambled_seqs[np.argmin(similarity)]
	return max_scrambled

def calculate_gc(seq):
	count = 0
	for i in range(len(seq)):
		if seq[i] == 'G' or seq[i] == 'C' or seq[i] == 'g' or seq[i] == 'c':
			count += 1

	return count / float(len(seq))


#------------------------------------------------------------------------
if __name__ == '__main__':
    
    # read in primers 5' to 3' of oligo
    fwd_primers = csv_reader('forward_primers.csv')
    rev_primers = csv_reader('reverse_primers.csv')

	# dictionary that will hold dictionaries of the subpool libraries
	# key will be name of subpool, easy to match with primers
    subpool_libraries = {}
    
    #import old backgrounds to be used in some designs
    old_bgs = csv_reader('oldbackgrounds.csv')

	# putting these in lowercase makes it easier to double check work
	# will convert to all uppercase for final version to be sent to Agilent
    consensus = 'TGACGTCA'.lower()
    consensus_flank = 'ATTGACGTCAGC'.lower()
    consensus_flank_rc = 'GCTGACGTCAAT'.lower()
    moderate = 'TGACGTCT'.lower()
    moderate_flank = 'ATTGACGTCTGC'.lower()
    weak = 'TGAAGTCA'.lower()
    weak_flank = 'ATTGAAGTCAGC'.lower()
    nosite = 'XXXXXXXX'.lower()
    nosite_flank = 'XXXXXXXXXXXX'.lower()


#------------------------------------------------------------------------------
#Generating the gc content background library with CREs

bgs_top3 = capitalize_values(csv_reader('bgs_top3.csv'))


#Add name for backgrounds with number and GC content

def name_background(backgrounds):
    bg_names = []
    name_number = 1
    seq = list(backgrounds.keys())
    seq_name = {}
    
    for bg in backgrounds.values():
        bg_name = 'background_' + str(name_number)
        name_number = name_number + 1
        bg_names.append(bg_name)
        
    for i in range(len(seq)):
        seq_name[bg_names[i]] = seq[i]
        
    return(seq_name)
    
bgs_top3_name = name_background(bgs_top3)

csv_writer(bgs_top3_name, 'bgs_top3_name.csv')


#Scramble 10 bp of the backgrounds in 5 bp increments to check for TFBSs 
#generated

tiles = {}

for x in bgs_top3_name:

	bg = bgs_top3_name[x]
	for i in range(0, len(bg), 5):
		if i + 10 > len(bg):
			continue
		to_scramble = bg[i:i+10]
		scrambled = most_scramble(to_scramble, 100)
		tile = bg[:i] + ''.join(scrambled).lower() + bg[i+10:]
		tile_name = 'scramble_dist_' + str(140-i) + '_similarity_' + str(similar(to_scramble, scrambled)) + '_' + x
		tiles[tile_name] = tile

	subpool_libraries['minP'] = tiles
	print(len(tiles), 'scrambled sequences in minP subpool')
    
csv_writer(tiles, 'cre_test_scramble.csv')


#------------------------------------------------------------------------------
#Add 0-6 consensus CRE sites to backgrounds

cre_vars = [nosite_flank, consensus_flank]
names = ['nosite', 'consensus']
cre_combs = {}
bg = []

for x in bgs_top3_name:
    bg = bgs_top3_name[x]
    # the six sites are in a constant position on the template, so parse out
    # the parts of the template we're going to use now
    spacers = [bg[i:i+13] for i in range(0, len(bg), 25)]
    
    for i in range(len(cre_vars)):
        for j in range(len(cre_vars)):
            for k in range(len(cre_vars)):
                for l in range(len(cre_vars)):
                    for m in range(len(cre_vars)):
                        for n in range(len(cre_vars)):
                            seq = spacers[0] + cre_vars[i] + spacers[1] + cre_vars[j]
                            seq += spacers[2] + cre_vars[k] + spacers[3] + cre_vars[l]
                            seq += spacers[4] + cre_vars[m] + spacers[5] + cre_vars[n] 
                            
                            # replace all  xxx no site with bases from template
                            loc = seq.find('xxxxxxxxxxxx')
                            while loc > -1:
                                seq = seq[:loc] + bg[loc:loc+12] + seq[loc+12:]
                                loc = seq.find('xxxxxxxxxxxx')
                                
                            gc = calculate_gc(seq)
                            
                            header = 'sixsite_' + '_'.join([names[i], names[j], 
                                       names[k], names[l], names[m], names[n], 
                                       x, 'gc', str(gc)])
        
                            cre_combs[header] = seq
                            subpool_libraries['minP'] = cre_combs
                            
        print(len(cre_combs), "affinity x background sequences in minP subpool")
        

csv_writer(cre_combs, 'cre_test_sixsite.csv')


#------------------------------------------------------------------------------
#Design the 2-site library again with CREs spaced 5-20 bp apart and CRE sites 
#varying between consensus CREs and their reverse complement. Move both CREs
#along backgrounds maintaining spacing. Keep in mind both CRE sites are 12 bp
#long

cre_spac_dist_orient = {}
bg = []
cre_sites = [consensus_flank, consensus_flank_rc]
cre_names = ['consensus', 'RCconsensus']

#spacing refers to number of bases in between the two CREs with flanks so 
#spacing + 4 refers to the actual number of bases between CREs

spacing = [1, 6, 11, 16, 66]

for x in bgs_top3_name:
    bg = bgs_top3_name[x]
    for space in spacing:
        for j in range(len(cre_sites)):
            for k in range(len(cre_sites)):
                seqs = [bg[:i] + cre_sites[j] + bg[i + len(cre_sites[j]):i + len(cre_sites[j]) + space] + cre_sites[k] + bg[i + 2*len(cre_sites[k]) + space:]
						for i in range(len(bg)-(2*12 + space) + 1)]

			twobs_spacing = {'2BS ' + str(space) + ' bp spacing consensus+flank x2_dist_' + str((150 - (2*len(consensus_flank)+space)-i)) + '_' + x : seqs[i]
								for i in range(len(seqs))}
			consensus_flank_spacing.update(twobs_spacing)

	subpool_libraries['subpool_3'] = consensus_flank_spacing
	print len(consensus_flank_spacing), "sequences in subpool 3"






#oldcode

#------------------------------------------------------------------------------
#Design an optimized library that tests 5 sites with spacing that places CREs
#either at one helical turn upstream or 2 helical turns. This results in 
#either alternating 2-3 bp spacing between the 8 bp CRE sites or 9 bp spacing
#between 12 bp CRE sites with flanks. From 5' to 3' site 1 and site 5 will just
#be consensus, 2 and 3 will be 1 of the 4 affinities and 4 will be either 
#moderate or consensus. Then each arrangement will then be moved in 1 bp
#increments to the 5' end of the oligo for 10 bp total movement. This will be
#performed on a subset of 3 backgrounds. 

site_2_3 = [nosite, weak, moderate, consensus]
site_2_3_names = ['nosite', 'weak', 'moderate', 'consensus']

site_4 = [moderate, consensus]
site_4_names = ['moderate', 'consensus']

spacing_2bp = {}

for x in optim_enhancers:
    enh = optim_enhancers[x]
    for j in range(len(site_2_3)):
        for k in range(len(site_2_3)):
            for l in range(len(site_4)):
                
                site1 = consensus
                site2 = site_2_3[j]
                site3 = site_2_3[k]
                site4 = site_4[l]
                site5 = consensus
                
                seqs = [enh[:i] + 'AT' + site1 + 'GCT' + site2 + 'GT' + site3 
                        + 'GAT' + site4 + 'GT' + site5 + 'GC' + enh[i+54:]
                for i in range(len(enh)- 54 - 10 + 1, len(enh)- 54 + 1, 1)]
                
                #loc = seqs.find('xxxxxxxx')
                #while loc > -1:
                    #seqs = seqs[:loc] + seqs[loc:loc+8] + seqs[loc+8]
                    #loc = seq.find('xxxxxxxx')
                
                spacing_name = {'optim_turn_1_consensus_' + names[j] + '_' +
                                names[k] + '_' + names[l] + '_consensus_dist_' 
                                + str(11 - i) + '_' + x : seqs[i] 
                                for i in range(len(seqs))}
                
                spacing_2bp.update(spacing_name)
                subpool_libraries['minP'] = spacing_2bp
                subpool_libraries['shortminP'] = spacing_2bp
                print(len(spacing_2bp), "optimized 1 turn sequences in both subpools")
    
csv_writer(spacing_2bp, 'cre_test_optim_1.csv')









def name_gc_content_background(backgrounds):
    gc_names = []
    name_number = 1
    seq = list(backgrounds.keys())
    seq_name = {}
    
    for gc in backgrounds.values():
        gc_name = 'background_' + str(name_number) + '_gc_' + str(gc)
        name_number = name_number + 1
        gc_names.append(gc_name)
        
    for i in range(len(seq)):
        seq_name[gc_names[i]] = seq[i]
        
    return(seq_name)
    
bgs_top4_name = name_gc_content_background(bgs_top4)


#This design will be present in both subpools to additionally test the effect 
#of distance to a shorter minimal promoter

















