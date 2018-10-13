# -*- coding: utf-8 -*-
"""
Functions for Johns Hopkins Coursera Bioinformatics Algorithms Course
GLobal Alignment Functions
Callen Hyland, November 2017
"""

from itertools import permutations


def edit_distance(x,y):
    """Create matrix of edit distances between two DNA sequences, x and y
    """
    D = []
    for i in range(len(x)+1):
        D.append([0]* (len(y)+1))
    
    # initialize first row and first column
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = i
    # go through rows and columns
    for i in range(1,len(x)+1):
        for j in range(1, len(y)+1):
            hdist = D[i][j-1]+1
            vdist = D[i-1][j]+1
            if x[i-1] == y[j-1]:
                ddist = D[i-1][j-1]
            else:
                ddist = D[i-1][j-1]+1
            D[i][j] = min(hdist,vdist,ddist)
    return D


def global_alignment(x,y, alphabet, score):
    """GLobally align two sequences with dynamic programming
    Input
        x : first sequence
        y : second sequence
        alphabet : DNA bases corresponding to columns and rows in score matrix
        score : score matrix for nucleotide substitutions
    Returns
        Matrix of scores for global alignment of sequences
    """
    D = []
    for i in range(len(x)+1):
        D.append([0]* (len(y)+1))
    
    # initialize first row and first column
    for i in range(1, len(x)+1):
        D[i][0] = D[0][i-1] + score[alphabet.index(x[i-1])][-1]
    for i in range(1, len(y)+1):
        D[0][i] = D[i-1][0] + score[-1][alphabet.index(y[i-1])]
        
    # go through rows and columns
    for i in range(1,len(x)+1):
        for j in range(1, len(y)+1):
            hdist = D[i][j-1] + score[-1][alphabet.index(y[j-1])]
            vdist = D[i-1][j] + score[alphabet.index(x[i-1])][-1]
            if x[i-1] == y[j-1]:
                ddist = D[i-1][j-1]
            else:
                ddist = D[i-1][j-1] + score[alphabet.index(x[i-1])][alphabet.index(y[j-1])]
            D[i][j] = min(hdist,vdist,ddist)
    return D[-1][-1]


def overlap(a,b,min_length = 3):
    """find the length of the longest overlap between a and b
    """
    start = 0
    while True:
        # find the first occurance of the prefix of b in a
        start = a.find(b[:min_length], start)
        # If the prefix of b does not occur anywhere in a, return zero overlap
        if start == -1:
            return 0
        # If there is some overlap, check that the rest of the suffix of a matches b beyond the min length
        if b.startswith(a[start:]):
            return len(a)-start
        # if the suffix of a is not a perfect match for the prefix of a, increment start and search again
        start += 1


def naive_overlap_map(reads,k):
    """
    Inputs
        reads : DNA sequences to be aligned
        k : minimum length of alignment
    Returns
        olaps : overlap map as dictionary
    """
    olaps = {}
    #loop through permuations of two read alignments
    for a,b in permutations(reads,2):
        #call overlap function to find the overlap between a and b
        olen = overlap(a, b, min_length = k)
        #if not zero, add length to dictionary
        if olen > 0:
            olaps[(a,b)] = olen
    return olaps


"""Script"""

# Standard DNA alphabet
alphabet = ['C','A','T','G']

# score matrix
score = [[0,4,2,4,8],\
         [4,0,4,2,8],\
         [2,4,0,4,8],\
         [4,2,4,0,8],\
         [8,8,8,8,8]]

x = 'ATCTACACTCGATGC'
y = 'TGCTACACCGATGC'

print(edit_distance(x,y))
print(global_alignment(x,y))
print(overlap(x,y))

