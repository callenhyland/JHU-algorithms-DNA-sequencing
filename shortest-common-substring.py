# -*- coding: utf-8 -*-
"""
Functions for Johns Hopkins Coursera Bioinformatics Algorithms Course
Shortest Common Substring
Callen Hyland, November 2017
"""

from itertools import permutations

def overlap(a,b,min_length = 3):
    #find the length of the longest overlap between a and b
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


def scs(ss):
    shortest_sup = None
    #loop through all permuations (orders) of substrings
    for ssperm in permutations(ss):
        #initialize superstring with the first string in the list
        sup = ssperm[0]
        #loop through strings in current order starting at zero
        for i in range(len(ss)-1):
            #find length of overlap between adjacent strings
            olen = overlap(ssperm[i],ssperm[i+1], min_length=1)
            #add the suffix of second string (after the overlap) onto the end of the superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup)<len(shortest_sup):
            shortest_sup = sup
            
    return shortest_sup


def sup_count(ss):
    shortest_sup = None
    sup_list = []
    #loop through all permuations (orders) of substrings
    for ssperm in permutations(ss):
        #initialize superstring with the first string in the list
        sup = ssperm[0]
        #loop through strings in current order starting at zero
        for i in range(len(ss)-1):
            #find length of overlap between adjacent strings
            olen = overlap(ssperm[i],ssperm[i+1], min_length=1)
            #add the suffix of second string (after the overlap) onto the end of the superstring
            sup += ssperm[i+1][olen:]
        sup_list.append(sup)
        if shortest_sup is None or len(sup)<len(shortest_sup):
            shortest_sup = sup

    sup_count = 0
    slen = len(shortest_sup)
    for sup in sup_list:
        if len(sup) == slen:
            sup_count += 1
            
    return sup_count


def pick_maximal_overlap(reads,k):
    #initialize reada, readb, and best_olen
    reada, readb = None, None
    best_olen = 0
    # loop through all pairs of two reads
    for a,b in permutations(reads,2):
        #find length of overlap
        olen = overlap(a,b,min_length=k)
        #is overlap better than current best overlap?
        if olen > best_olen:
            reada = a
            readb = b
            best_olen = olen
            
    return reada, readb, best_olen


def greedy_scs(reads, k):
    #initialie by picking the first maximal overlap
    reada, readb, olen = pick_maximal_overlap(reads,k)
    while olen > 0:
        #remove reada and readb and replace with the result of their overlap
        reads.remove(reada)
        reads.remove(readb)
        #concatenate reada with the suffix of readb
        reads.append(reada + readb[olen:])
        #find the maximal overlap in the new set
        reada, readb, olen = pick_maximal_overlap(reads,k)
    #concatenate all remaining reads that do not have overlap
    return ''.join(reads)


def de_bruijnize(st,k):
    #initialize list of edges and set of nodes
    edges = []
    nodes = set()
    #loop through all kmers in the string
    for i in range(len(st) - k +1):
        #add the pair of k-1mers as a tuple to the list of edges
        edges.append(st[i:i+k-1],st[i+1:i+k])
        #add the nodes to the set of nodes
        nodes |= set(st[i+1:i+k])
        nodes |= set(st[i:i+k-1])
    return edges, nodes


def overlap_map_scs(reads,k):
    #initialize empty dictionary to hold sets of reads matching kmers
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
    
    olaps = []
    olens = []
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
                    olaps.append((read,m))
                    olens.append(olen)
                    #edge = True

    return olaps, olens    #, edge_count


def fast_greedy_scs(reads, k):
    #initialie by picking the first maximal overlap
    #reada, readb, olen = pick_maximal_overlap(reads,k)
    olaps, olens = overlap_map_scs(reads,k)
    max_olen = 0
    for i in range(len(olens)):
        if olens[i] > max_olen:
            max_olen = olens[i]
            reada = olaps[i][0]
            readb = olaps[i][1]
    
    j=0
    while max_olen > 0:
        j+=1
        #print(j)
        #remove reada and readb and replace with the result of their overlap
        reads.remove(reada)
        reads.remove(readb)
        #concatenate reada with the suffix of readb
        reads.append(reada + readb[max_olen:])
        #print(reada + readb[max_olen:])
        #find the maximal overlap in the new set
        #reada, readb, olen = pick_maximal_overlap(reads,k)
        olaps, olens = overlap_map_scs(reads,k)
        max_olen = 0
        for i in range(len(olens)):
            if olens[i] > max_olen:
                max_olen = olens[i]
                reada = olaps[i][0]
                readb = olaps[i][1]
    #concatenate all remaining reads that do not have overlap
    return ''.join(reads)