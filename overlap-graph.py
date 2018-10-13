# -*- coding: utf-8 -*-
"""
Functions for Johns Hopkins Coursera Bioinformatics Algorithms Course
Boyer-Moore algorithm
Callen Hyland, November 2017
"""

def overlap(a,b,min_length = 3):
    """Finds the length of the longest overlap between a and b
    """
    start = 0
    while True:
        #find the first occurance of the prefix of b in a
        start = a.find(b[:min_length], start)
        #If the prefix of b does not occur anywhere in a, return zero overlap
        if start == -1:
            return 0
        #If there is some overlap, check that the rest of the suffix of a matches b beyond the min length
        if b.startswith(a[start:]):
            return len(a)-start
        #if the suffix of a is not a perfect match for the prefix of a, increment start and search again
        start += 1


def overlap_map(reads,k):
    """Finds all overlap graph for a list of reads
    Inputs
        reads : list of DNA sequences
        k : length of kmers to split the reads up into
    Returns
        olaps : overlap graph as dictionary
        edge_count : number of edges in overlap graph
    """
    dt = {}
    #loop through all of the reads and 
    for read in reads:
        #find all kmers in read
        for i in range(len(read)-k+1):
            kmer = read[i:i+k]
            #if the kmer is already in the dictionary, add the read to the set of matches
            if kmer in dt:
                s = dt[kmer]
                s.add(read)
                dt[kmer] = s
            #if the kmer is not already in the dictionary, add it as a key and create a set containing the read
            else:
                dt[kmer] = {read}
    
    olaps = {}
    edge_count = 0
    
    #loop through all reads again
    for read in reads:
        edge = False
        #extract k-length suffix from each read
        suf = read[len(read)-k:]
        #find the set of reads corresponding that contain that kmer
        matches = dt[suf]
        #loop through reads in set
        for m in matches:
            #check that we are not matching the read with itself
            if m != read:
                olen = overlap(read,m,min_length=k)
                if olen > 0:
                    olaps[(read,m)] = olen
                    edge = True
        if edge:
            edge_count += 1

    return olaps, edge_count