# -*- coding: utf-8 -*-
"""
Functions for Johns Hopkins Coursera Bioinformatics Algorithms Course
Boyer-Moore algorithm
Callen Hyland, November 2017
"""

from bm_preproc import BoyerMoore

def boyer_moore_with_counts(p, p_bm, t): 
    """ Do Boyer-Moore matching with counts extending boyer_moore function
    Inputs
        p : pattern to match against dna sequence
        t : dna sequence
        p_bm : BoyerMoore object for p 
    Returns:
        occurences : numper of time p appears in t
        num_alignments : number of alignments performed
        num_character_comparisons: total number of characters compared
    """
    i = 0
    occurrences = []
    num_alignments = 0
    num_character_comparisons = 0
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        num_alignments += 1
        for j in range(len(p)-1, -1, -1):
            num_character_comparisons += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return(occurrences, num_alignments, num_character_comparisons)


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


"""Script to run boyer_moore_algorithm on human chromosome and given pattern"""

p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
t = readGenome('chr1.GRCh38.excerpt.fasta')
p_bm = BoyerMoore(p, 'ATGC')
matches, alignments, char_comparisons = boyer_moore_with_counts(p, p_bm, t)
print(matches, alignments, char_comparisons)

