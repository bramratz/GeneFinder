#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 07:25:04 2020

@author: bramratz
"""

# ---------------------------------------------------------------------

### Dependencies 

# ---------------------------------------------------------------------

import random

# ---------------------------------------------------------------------

### Loading data

# ---------------------------------------------------------------------

def loadSeq(filename):
    """
    Given a fasta file returns the sequence associated
    with the ID
    """
    seqList = [] # Create list to save sequence to
    with open(filename) as f:
        # Save lines in fasta as list 
        lines = [item.strip() for item in f.readlines() if not item == '']
        for line in lines:
            if line.startswith('>'):
                pass # '>' are ID lines, so skip
            else:
                seqList.append(line) # Append all other lines as these represent sequences to main list 
    return "".join(seqList) # Return list joined as a string 

### Base Compliments ###
def compBase(N: str) -> str:
    """
    Given a nucleotide base (A, C, G, T) returns the 
    compliment of that base
    
    Given A -> Returns T
    """
    # Dictionary of complimentary bases
    baseDict = {'A': 'T', 
                'C': 'G',
                'G': 'C',
                'T': 'A'}
    return baseDict[N] # Return compliment

### Reverse Compliments ###
def reverse(s: str) -> str:
    """
    Given a DNA sequence, s, returns a string that is
    the reverse of s
    
    Given "AATT" -> Returns "TTAA"
    """
    compliment = [] # Initialize list to hold bases in reverse order
    # Iterate string in reverse, end to start
    for i in range(len(s)-1, -1, -1):
        compliment.append(s[i]) # Add base to compliment list 
    return "".join(compliment) # return compliment 

def reverseCompliment(DNA: str) -> str:
    """
    Given a DNA sequence (5' -> 3') returns the reverse
    compliment. ie. The compliment of the DNA seq in 
    (5' -> 3')
    
    Given "TTGAC" ->  Returns "GTCAA"
    """
    compliment = '' # Initialize empty string to hold complimentary seq
    for i in DNA:
        compliment += compBase(i) # Append complimentary base to string
    
    return reverse(compliment) # Return the reverse of the complimentary string


### Finding ORFs ###

def restOfORF(DNA: str) -> str:
    """
    Given a DNA sequence finds the first in frame stop
    codon and returns the sequence from the start up 
    to, but not including the stop codon. This is an 
    ORF.
    Assumes that the DNA sequence is 5'->3' and 
    that it begins aith a start codon, 'ATG'.
    If no stop codon found returns the whole sequence
    """
    ORF = '' # empty string for ORF
    for i in range(0, len(DNA)-2, 3):
        # if codon isn't a stop codon append the codon to ORF
        if DNA[i:i+3] not in ['TAG', 'TAA', 'TGA']:
            ORF += DNA[i:i+3]
        # Stop if a stop codon is found
        elif DNA[i:i+3] in ['TAG', 'TAA', 'TGA']:
            break 
    return ORF

def oneFrameV2(DNA: str) -> str:
    """
    Given a DNA sequence searches from position 0, 3
    bases at a time for a start codon, 'ATG'. Calls
    restOfORF to find the ORF associated with that start
    codon and stores the sequence in a list.
    Ignores nested start codons. If nested start codons
    are present will only return the longest one.
    Returns a list of all ORF found in the DNA sequence.
    """
    frames = [] # empty list for ORFs found
    frameIdx = 0  # length of DNA seq to iterate
    while frameIdx < len(DNA) - 2:
        if DNA[frameIdx:frameIdx+3] == "ATG":
            ORF = restOfORF(DNA[frameIdx:])
            frames.append(ORF)
            frameIdx += len(ORF)
        else:
            frameIdx += 3
    return frames 

def longestORFV2(DNA: str) -> str:
    """
    Given a DNA sequence returns the sequence of
    the longest ORF in any of the possible 3 
    reading frames.
    """
    # Dictionary to hold ORF's in each frame
    frames = {'frame1': '',
              'frame2': '',
              'frame3': ''}
    # Find ORF's for three possible reading frames
    frames['frame1'] = oneFrameV2(DNA)
    frames['frame2'] = oneFrameV2(DNA[1:])
    frames['frame3'] = oneFrameV2(DNA[2:])
    
    length = 0 # length of sequence
    longestSeq = '' # DNA Sequence
    for key, val in frames.items():
        for seq in val:
            # Update longestSeq only if sequence is
            #   longer than previous
            if len(seq) > length:
                length = len(seq)
                longestSeq = seq
    return longestSeq

### Find longest ORF on either strand of DNA ###

def longestORFBothStrands(DNA: str) -> str:
    """
    Given a DNA sequence finds the longest ORF in both the 
    5' -> 3' and the 3' -> 5' directions. Returns the
    longest ORF (tie breaks arbitrarily).
    """
    regStrand = longestORFV2(DNA) # Find longest ORF of DNA seq
    revDNA = reverseCompliment(DNA) # Get the reverse compliment of the DNA strand
    revComp = longestORFV2(revDNA) # Find longest ORD of the reverse compliment 
    
    if len(regStrand) > len(revComp):
        return regStrand # Return ORF of DNA seq if longer than reverse compliment 
    elif len(regStrand) == len(revComp):
        return regStrand # Return ORF of DNA seq same length as reverse compliment
    else:
        return revComp # Return ORF of reverse compliment seq if longer than DNA seq


### Assess whether ORF are genes or likely to appear randomly ###

def longestORFNoncoding(DNA: str, numReps: int) -> int:
    """
    Given a DNA sequence shuffles the sequence to create
    a random sequence and finds the longest ORF. Repeats 
    this process numReps number of times. The longest 
    ORF of all the iterations is returned. Idea is to 
    compare to the longest ORF found in the original 
    sequence to determine is the sequence length 
    indicative of a coding region or would I expect to see
    sequences this long in random, non-coding sequences
    """
    longestORF = 0 # Holds length of longest ORF found
    reps = 0 # Counts iterations
    
    # Iterate numReps times 
    while numReps > reps:
        # Shuffle the sequence
        seqList = list(DNA) # Convert str to list 
        random.shuffle(seqList) # shuffle the list in place
        shuffledSeq = ''.join(seqList) # Convert back to str
    
        # Calculate the longest ORF
        shuffledSeqORF = longestORFBothStrands(shuffledSeq)
        
        # Compare to current longest
        if len(shuffledSeqORF) > longestORF:
            longestORF = len(shuffledSeqORF)
        
        # Increase counter
        reps += 1
    
    return longestORF

### ID real ORFs ###

def findORFs(DNA: str):
    """
    Given a DNA sequence identifies all real ORFs in the 
    sequence in the 5' -> 3' directions for each reading
    frame and returns them as a list, if there are none
    returns an empty list
    """
    # Find ORF's for three possible reading frames
    frame1 = oneFrameV2(DNA)
    frame2 = oneFrameV2(DNA[1:])
    frame3 = oneFrameV2(DNA[2:])
    
    # Combine and return as a list 
    return frame1 + frame2 + frame3

def findORFsBothStrands(DNA: str):
    """
    Similar to findORFs, given a DNA sequence identifies
    all ORFs in the sequence, but in both the forward 
    (5'->3') and reverse (3'->5') driections. Returns a 
    list of all the ORFs found 
    """
    revDNA = reverseCompliment(DNA) # Get the reverse compliment of the DNA strand
    
    # Find ORF's for three possible reading frames in forward direction
    frame1 = oneFrameV2(DNA)
    frame2 = oneFrameV2(DNA[1:])
    frame3 = oneFrameV2(DNA[2:])
    
    # Repeat for the reverse compliment 
    revFrame1 = oneFrameV2(revDNA)
    revFrame2 = oneFrameV2(revDNA[1:])
    revFrame3 = oneFrameV2(revDNA[2:])
    
    # Combine and return as a list 
    return frame1 + frame2 + frame3 + revFrame1 + revFrame2 + revFrame3

### Finding coordinates ###

def getCoordinates(orf: str, DNA: str):
    """
    For a given ORF returns the coordinates [start, end]
    of the ORF in the DNA sequence
    """
    # Find coordinates
    start = DNA.find(orf)
    
    if start == -1:
        # Not present on the forward strand. Try reverse strand
        revORF = reverseCompliment(orf) # Get the reverse compliment of the orf
        start = DNA.find(revORF)
        if start == -1:
            # Not present on reverse compliment either
            return "ORF not present in given DNA strand"
        
        return [start, start + len(orf)] # Present on reverse strand
    return [start, start + len(orf)] # Present on forward strand

### Gene finding ###

def geneFinder(DNA: str, minLen: int):
    """
    Given a DNA sequence identifies ORFs longer than minLen
    returning list where each item is a list containing
    [startCorrdinate, endCoordinate, proteinSeq]. The 
    list is sorted from smalledst starting to coordinate
    to largest. ie. in order found in sequence.
    """
    # Find all ORFs on both strands
    allORFs = findORFsBothStrands(DNA)
    
    # Filter by size and calculate coordinates + aa seq
    finalORF = [] # List to hold information about ORFs that are larger than the minLen
    for i in allORFs:
        if len(i) > minLen:
            coord = getCoordinates(i, DNA) # Get ORF coordinates
            aa = codingStrandToAA(i) # Get ORF amino acid sequence
            
            # Append results to final list 
            finalORF.append([coord[0], coord[1], aa])
            
    finalORF.sort() # Order the list by start coordinate
    return finalORF