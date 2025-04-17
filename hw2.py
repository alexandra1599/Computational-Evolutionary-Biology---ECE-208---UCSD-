#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 14:46:32 2023

@author: alexandramikhael
"""

from math import log2
import csv
from threeway_align.BLOSUM62 import *
AMINOS   = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

def compute_gap_penalty(indel_rate,P):
    '''
    This function computes the gap penalty for threeway_align using the BLOSUM62 model and the indel rate
    The transition probabilities have been loaded for you and stored in the nested dictionary P (see file threeway_align/BLOSUM62 and threeway_align/blosum62sym.csv ) 
    '''
    gap = int(2*math.log2(indel_rate*19/20))
    # TODO: replace with your way to compute gap penalty
    return gap

def threeway_align(s1, s2, s3, B, gap, VERBOSE=False):
    '''
    This function computes the threeway sequence alignment between strings s1, s2, and s3 given a substitution matrix and gap penalty
    :param s1 (string): the first sequence
    :param s2 (string): the second sequence
    :param s3 (string): the third sequence
    :param B (char,char -> int) : the substitution matrix (e.g. BLOSUM62)
    :param gap (negative float):  gap penalty
    '''
    
    sz_s1 = len(s1)+1
    sz_s2 = len(s2)+1
    sz_s3 = len(s3)+1
    # initialize (S[i][j][k] = (score,arrow))
    if VERBOSE:
        from sys import stderr; print("Initializing cube", file=stderr)
    S = [[[None for k in range(len(s3)+1)] for j in range(len(s2)+1)] for i in range(len(s1)+1)]
    S[0][0][0] = (0, None)        # corner

    # fill in cube axes
    if VERBOSE:
        print("Filling cube axes", file=stderr)
    for i in range(1, len(s1)+1): # s1 axis
        S[i][0][0] = (i,None); # TODO: replace with your code
    for j in range(1, len(s2)+1): # s2 axis
        S[0][j][0] = (j,None); # TODO: replace with your code
    for k in range(1, len(s3)+1): # s3 axis
        S[0][0][k] = (k,None); # TODO: replace with your code

    # fill in cube faces
    if VERBOSE:
        print("Filling cube faces", file=stderr)
    # TODO: add your code here
    #filling the k=0 face
    for i in range(1,len(s1)+1):
        for j in range(1,len(s2)+1):
           S[i][j][0] = max(
                    S[i-1][j-1][0] + B[AMINOS.index(s1[i-1])][AMINOS.index(s2[j-1])]+gap,
                    S[i-1][j][0] + gap,
                    S[i][j-1][0] + gap
                    )
    #filling the j=0 face
    for i in range(1,len(s1)+1):
        for k in range(1,len(s3)+1):
             S[i][0][k] = max(
                    S[i-1][0][k-1] + B[AMINOS.index(s1[i-1])][AMINOS.index(s3[k-1])]+gap,
                    S[i-1][0][k] + gap,
                    S[i][0][k-1] + gap
                    )
    #filling the i=0 face
    for j in range(1,len(s2)+1):
        for k in range(1,len(s3)+1):
             S[0][j][k] = max(
                    S[0][j-1][k-1] + B[AMINOS.index(s1[i-1])][AMINOS.index(s3[k-1])]+gap,
                    S[][j-1][k] + gap,
                    S[0][j][k-1] + gap
             )

    # fill in rest of cube
    if VERBOSE:
        print("Filling rest of cube", file=stderr)
    # TODO: add your code here
    diff = 2*max( #set a maximum amount of deviation from the diagonal in any direction
            abs(sz_s1-sz_s2),
            abs(sz_s2-sz_s1),
            abs(sz_s3-sz_s2),
            5
            )
    for i in range (1,sz_s1):
        #print(i)
        for j in range (1,sz_s2):
            if (abs(i-j)< diff):
                for k in range (1,sz_s3):
                    if (abs(i-k)< diff) and (abs(j-k)< diff):
                        S[i,j,k] = max(
                            S[i,j-1,k] + gap, #1
                            S[i,j-1, k-1] + B[AMINOS.index(s2[j-1])][AMINOS.index(s3[k-1])]+ gap, #2
                            S[i, j, k-1] + gap, #3
                            S[i-1, j-1, k] + B[AMINOS.index(s1[i-1])][AMINOS.index(s2[j-1])] + gap, #4
                            S[i-1,j-1,k-1] + B[AMINOS.index(s1[i-1])][AMINOS.index(s2[j-1])] + B[AMINOS.index(s2[j-1])][AMINOS.index(s3[k-1])] + B[AMINOS.index(s1[i-1])][AMINOS.index(s3[k-1])], #5
                            S[i-1,j, k-1] + B[AMINOS.index(s1[i-1])][AMINOS.index(s3[k-1])] + gap,#6
                            S[i-1, j, k] + gap, #7
                            )

    # backtrack to get alignments
    x = sz_s1-1
    y = sz_s2-1
    z = sz_s3-1
    if VERBOSE:
        print("Backtracking to build alignment", file=stderr)
    aln_s1 = ""; aln_s2 = ""; aln_s3 = ""
    
    #check all 7 movement posibilities and figure out which one happened
    while x > 0 and y > 0 and z > 0:
        xyzCell = S[x,y,z]
        #print(x,y,z)
        if xyzCell == S[x,y-1,z] + gap: #1
            aln_s1 = '-'+aln_s1
            aln_s2 = s2[y-1]+aln_s2
            aln_s3 = '-'+aln_s3
            y=y-1
        elif xyzCell == S[x,y-1, z-1] + B[AMINOS.index(s2[y-1])][AMINOS.index(s3[z-1])]+ gap: #2
            aln_s1 = '-'+aln_s1
            aln_s2 = s2[y-1]+aln_s2
            aln_s3 = s3[z-1]+aln_s3
            y=y-1
            z=z-1
        elif xyzCell == S[x, y, z-1] + gap: #3
            aln_s1 = '-'+aln_s1
            aln_s2 = '-'+aln_s2
            aln_s3 = s3[z-1]+aln_s3
            z=z-1
        elif xyzCell == S[x-1, y-1, z] + B[AMINOS.index(s1[x-1])][AMINOS.index(s2[y-1])] + gap: #4
            aln_s1 = s1[x-1]+aln_s1
            aln_s2 = s2[y-1]+aln_s2
            aln_s3 = '-'+aln_s3
            y=y-1
            x=x-1
        elif xyzCell == S[x-1,y-1,z-1] + B[AMINOS.index(s1[x-1])][AMINOS.index(s2[y-1])] + B[AMINOS.index(s2[y-1])][AMINOS.index(s3[z-1])] + B[AMINOS.index(s1[x-1])][AMINOS.index(s3[z-1])]: #5
            aln_s1 = s1[x-1]+aln_s1
            aln_s2 = s2[y-1]+aln_s2
            aln_s3 = s3[z-1]+aln_s3
            y=y-1
            x=x-1
            z=z-1
        elif xyzCell == S[x-1,y, z-1] + B[AMINOS.index(s1[x-1])][AMINOS.index(s3[z-1])] + gap:#6
            aln_s1 = s1[x-1]+aln_s1
            aln_s2 = '-'+aln_s2
            aln_s3 = s3[z-1]+aln_s3
            x=x-1
            z=z-1
        elif xyzCell == S[x-1, y, z] + gap: #7
            aln_s1 = s1[x-1]+aln_s1
            aln_s2 = '-'+aln_s2
            aln_s3 = '-'+aln_s3
            x=x-1
            
    while x > 0 and y > 0:
        xyzCell = S[x,y,z]
        if xyzCell == S[x,y-1,z] + gap: #1
            aln_s1 = '-'+aln_s1
            aln_s2 = s2[y-1]+aln_s2
            aln_s3 = '-'+aln_s3
            y=y-1
        elif xyzCell == S[x-1, y-1, z] + B[AMINOS.index(s1[x-1])][AMINOS.index(s2[y-1])] + gap: #4
            aln_s1 = s1[x-1]+aln_s1
            aln_s2 = s2[y-1]+aln_s2
            aln_s3 = '-'+aln_s3
            y=y-1
            x=x-1
        elif xyzCell == S[x-1, y, z] + gap: #7
            aln_s1 = s1[x-1]+aln_s1
            aln_s2 = '-'+aln_s2
            aln_s3 = '-'+aln_s3
            x=x-1
    while x > 0 and z > 0:
        xyzCell = S[x,y,z]
        if xyzCell == S[x, y, z-1] + gap: #3
            aln_s1 = '-'+aln_s1
            aln_s2 = '-'+aln_s2
            aln_s3 = s3[z-1]+aln_s3
            z=z-1
        elif xyzCell == S[x-1,y, z-1] + B[AMINOS.index(s1[x-1])][AMINOS.index(s3[z-1])] + gap:#6
            aln_s1 = s1[x-1]+aln_s1
            aln_s2 = '-'+aln_s2
            aln_s3 = s3[z-1]+aln_s3
            x=x-1
            z=z-1
        elif xyzCell == S[x-1, y, z] + gap: #7
            aln_s1 = s1[x-1]+aln_s1
            aln_s2 = '-'+aln_s2
            aln_s3 = '-'+aln_s3
            x=x-1
    while y > 0 and y > 0:
        xyzCell = S[x,y,z]
        if xyzCell == S[x,y-1,z] + gap: #1
            aln_s1 = '-'+aln_s1
            aln_s2 = s2[y-1]+aln_s2
            aln_s3 = '-'+aln_s3
            y=y-1
        elif xyzCell == S[x,y-1, z-1] + B[AMINOS.index(s2[y-1])][AMINOS.index(s3[z-1])]+ gap: #2
            aln_s1 = '-'+aln_s1
            aln_s2 = s2[y-1]+aln_s2
            aln_s3 = s3[z-1]+aln_s3
            y=y-1
            z=z-1
        elif xyzCell == S[x, y, z-1] + gap: #3
            aln_s1 = '-'+aln_s1
            aln_s2 = '-'+aln_s2
            aln_s3 = s3[z-1]+aln_s3
            z=z-1
    while z > 0:
        aln_s1 = '-'+aln_s1
        aln_s2 = '-'+aln_s2
        aln_s3 = s3[z-1]+aln_s3
        z=z-1
    while y > 0:
        aln_s1 = '-'+aln_s1
        aln_s2 = s2[y-1]+aln_s2
        aln_s3 = '-'+aln_s3
        y=y-1
    while x > 0:
        aln_s1 = s1[x-1]+aln_s1
        aln_s2 = '-'+aln_s2
        aln_s3 = '-'+aln_s3
        x=x-1 
 
    score = S[-1][-1][-1]
    return aln_s1[::-1],aln_s2[::-1],aln_s3[::-1],score 

def threeway_align_indel_rate(s1, s2, s3,B,indel_rate,VERBOSE=False):
    gap = compute_gap_penalty(indel_rate,P)                            
    a1,a2,a3,score = threeway_align(s1, s2, s3,B,gap, VERBOSE=VERBOSE) # EXTRA CREDITS: replace with your own algorithm to do alignment
    return a1,a2,a3,score  # optional (extra credits): replace with your own way to do alignment