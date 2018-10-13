# -*- coding: utf-8 -*-
"""
Functions for Johns Hopkins Coursera Bioinformatics Algorithms Course
Callen Hyland, November 2017
"""

def read_fastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline() # skip name line
            seq = fh.readline().rstrip() # read base sequence
            fh.readline() # skip placeholder line
            qual = fh.readline().rstrip() #base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


def findGCByPos(reads):
    ''' Finds the GC ratio at each position in the read 
    '''
    gc = [0] * 100
    totals = [0] * 100
    for read in reads:
        for i in range(len(read)):
            if read[i] == 'C' or read[i] == 'G':
                gc[i] += 1
            totals[i] += 1
    # Divide G/C counts by total counts to get the average at each position
    for i in range(len(gc)):
        if totals[i] > 0:
            gc[i] /= float(totals[i])
    return gc

def naive(t, p):
    '''Simple function for finding how many positions are matched
    between two DNA sequences t and p
    '''
    matches = []
    for i in range(len(t)-len(p)+1):
        match = True
        for j in range(len(p)):
            if not t[i+j] == p[j]:
                match = False
                break
        if match:
            matches.append(i)   
    return matches


def naive_mismatch(t, p, mm):
    '''Simple function for finding how many positions are matched
    between two DNA sequences t and p, allows up to mm mismatches
    '''
    matches = []
    for i in range(len(t)-len(p)+1):
        match = True
        mis = 0
        for j in range(len(p)):
            if not t[i+j] == p[j]:
                mis += 1
            if mis > mm:
                match = False
                break
        if match:
            matches.append(i)
    return matches


def longest_common_prefix(s1, s2):
    i = 0
    while i < len(s1) and i < len(s2) and s1[i] == s2[i]:
        i+=1  
    return s1[:i]


def read_genome(filename):
    genome = ''
    with open(filename,'r') as f:
        for line in f:
            if not line[0]== '>':
                genome += line.rstrip()
    return genome


def phred33toQ(qual):
    return ord(qual)-33


def create_hist(qualityStrings):
    # Create a histogram of quality scores
    hist = [0]*50
    for read in qualityStrings:
        for phred in read:
            q = phred33toQ(phred)
            hist[q] += 1
    return hist

