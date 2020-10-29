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
import sys
from aminoAcids import *


# ---------------------------------------------------------------------

### Loading data

# ---------------------------------------------------------------------

# Loads sequences from FASTA file one sequence at a time
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

# -------------------------------------------------------------------
    
### Working with double stranded DNA
    
# -------------------------------------------------------------------

# Find the complimentary base of a given nucleotide
def compBase(N: str) -> str:
    """
    Given a nucleotide base (A, C, G, T) returns the compliment of that 
    base. ie:
    
    Given A -> Return T
    """
    # Dictionary of complimentary bases
    baseDict = {'A': 'T', 
                'C': 'G',
                'G': 'C',
                'T': 'A'}
    
    return baseDict[N] 

# Find the reverse compliment of a nucleotide string
def reverseCompliment(DNA: str) -> str:
    """
    Given a DNA sequence (5' -> 3') returns the reverse
    compliment of the DNA seq in the (5' -> 3') direction. ie:
    
    Given "TTGAC" ->  Returns "GTCAA"
    """
    compliment = '' # To hold complimentary seq
    
    for i in DNA:
        compliment += compBase(i) # Add complimentary base to string
    
    return compliment[::-1] # Return the seq in reverse

assert reverseCompliment('AACC') == 'GGTT'

# -----------------------------------------------------------------

# Converting DNA sequences to amino acid sequences

# -----------------------------------------------------------------

# Convert a DNA codon to an amino acid
def amino(codon: str) -> str:
    """
    Given a string of three bases (codon) returnsthe corresponding amino 
    acid.
    """
    # Only 20 aa, so iterate through 0 - 20
    for i in range(0, 21, 1):
        for val in codons[i]: # Iterate codon or codons that correspond to aa at that index
            if codon == val: # Stop if match found
                return aa[i]

# Convert a DNA sequence to an amino acid sequence
def codingStrandToAA(DNA: str) -> str:
    """
    Given a DNA sequence returns the corresponding string of amino acids
    """
    AAstring = ''
    
    # Iterate length of DNA sequence 3 at a time (codon at a time)
    for i in range(0, len(DNA)-2, 3):
        AAstring += amino(DNA[i:i+3]) # Call amino function to determine correct aa for codon
    
    return AAstring   

# ------------------------------------------------------------------
    
### Finding ORFs

# ------------------------------------------------------------------

# Finds stop codons. Returns ORFs between start codon and stop codon.
def restOfORF(DNA: str) -> str:
    """
    Given a DNA sequence finds the first in frame stop codon and returns 
    the sequence from the start up to, but not including the stop codon. 
    This is an ORF. Assumes that the DNA sequence is 5'->3' and 
    that it begins aith a start codon, 'ATG'. If no stop codon found 
    returns the whole sequence.
    """
    ORF = '' # To hold ORF str
    
    # Iterate with step set to 3 since 3 bases per codon
    for i in range(0, len(DNA)-2, 3):
        # if codon isn't a stop codon append the codon to ORF
        if DNA[i:i+3] not in ['TAG', 'TAA', 'TGA']:
            ORF += DNA[i:i+3]
        # Stop if a stop codon is found
        elif DNA[i:i+3] in ['TAG', 'TAA', 'TGA']:
            break 
    
    return ORF

# Find ORFs in a DNA sequence
def oneFrameV2(DNA: str) -> str:
    """
    Given a DNA sequence searches from position 0, 3 bases at a time for a 
    start codon, 'ATG'. Calls restOfORF to find the ORF associated with that 
    startcodon and stores the sequence in a list. Ignores nested start 
    codons. If nested start codons are present will only return the longest 
    one. Returns a list of all ORF found in the DNA sequence.
    """
    frames = [] # List for ORFs found
    frameIdx = 0  # Counter for string index
    
    while frameIdx < len(DNA) - 2:
        
        if DNA[frameIdx:frameIdx+3] == "ATG":
            ORF = restOfORF(DNA[frameIdx:]) # Call restOfORF if start codon found
            frames.append(ORF) # Appends the ORF to frames list
            frameIdx += len(ORF) # Increases index by length of ORF
        else:
            frameIdx += 3 # Increases index by codon length otherwise
    return frames 

def longestORFV2(DNA: str) -> str:
    """
    Given a DNA sequence returns the sequence of the longest ORF in any of the 
    possible 3 reading frames.
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
            
            if len(seq) > length:
                length = len(seq) 
                longestSeq = seq # Update only if sequence is longer than previous
    
    return longestSeq

# Find longest ORF on either strand of DNA
def longestORFBothStrands(DNA: str) -> str:
    """
    Given a DNA sequence finds the longest ORF in both the 5' -> 3' and the 
    3' -> 5' directions. Returns the longest ORF (tie breaks to the sense 
    strand).
    """
    regStrand = longestORFV2(DNA) # Find longest ORF of DNA seq
    revDNA = reverseCompliment(DNA) # Get reverse compliment of the DNA strand
    revComp = longestORFV2(revDNA) # Find longest ORF of the reverse compliment 
    
    if len(regStrand) > len(revComp):
        return regStrand # Return ORF of DNA seq if longer than reverse compliment 
    
    elif len(regStrand) == len(revComp):
        return regStrand # Return ORF of DNA seq if same length as reverse compliment
    
    else:
        return revComp # Return ORF of reverse compliment seq if longer than DNA seq

# ------------------------------------------------------------------
        
### Comparing found ORF to random sequences 
        
# -----------------------------------------------------------------

# Compares ORF to a random sequence. Determines if you'd expect an ORF at 
#   least that long by chance
def longestORFNoncoding(DNA: str, numReps: int) -> int:
    """
    Given a DNA sequence shuffles the sequence to createa random sequence and 
    finds the longest ORF. Repeats this process numReps number of times. The 
    longest ORF of all the iterations is returned. Idea is to compare to the 
    longest ORF found in the original sequence to determine is the sequence 
    length indicative of a coding region or would you expect to see sequences 
    this long in random, non-coding sequences?
    """
    longestORF = 0 # longest ORF found
    reps = 0 # Counter
    
    # Iterate numReps times 
    while numReps > reps:
        seqList = list(DNA) # Convert str to list 
        random.shuffle(seqList) # Shuffle the list in place for new order
        shuffledSeq = ''.join(seqList) # Convert back to str
    
        # Calculate the longest ORF in shuffled sequence
        shuffledSeqORF = longestORFBothStrands(shuffledSeq)
        
        # Compare to current longest
        if len(shuffledSeqORF) > longestORF:
            longestORF = len(shuffledSeqORF)
        
        reps += 1 # Increase counter
    
    return longestORF

# -------------------------------------------------------------
    
### Gene finding
    
# ------------------------------------------------------------

# Find ORF on both strands in all frames
def findORFsBothStrands(DNA: str):
    """
    Given a DNA sequence identifies all ORFs in the sequence, in both the 
    forward (5'->3') and reverse (3'->5') driections. Returns a list of all 
    the ORFs found. 
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

# Finding coordinates of an ORF in a sequence 
def getCoordinates(orf: str, DNA: str):
    """
    For a given ORF returns the coordinates [start, end] of the ORF in the 
    DNA sequence
    """
    # Find coordinates
    start = DNA.find(orf)
    
    if start == -1: # Not present on the forward strand. Try reverse strand
        revORF = reverseCompliment(orf) # Get the reverse compliment of the orf
        start = DNA.find(revORF) # Find ORFs on reverse compliment 
        
        if start == -1: # Not present on reverse compliment either
            return "ORF not present in given DNA strand"
        
        return [start, start + len(orf)] # Present on reverse strand
    return [start, start + len(orf)] # Present on forward strand

# Find True ORFs and return their location in the sequence
def geneFinder(DNA: str, minLen: int):
    """
    Given a DNA sequence identifies ORFs longer than minLen returning list 
    where each item is a list containing [startCorrdinate, endCoordinate, 
    proteinSeq]. The list is sorted from smalledst starting to coordinate
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

# ----------------------------------------------------------
    
### Main
    
# ----------------------------------------------------------

IDList = [] # List of sequence IDs
seqList = [] # List of sequences

# Open and parse input file
with open(sys.argv[1]) as f:
    lines = [item.strip() for item in f.readlines() if not item == '']
    tempList = [] # List for sequences associated with each ID. Solves sequences on multiple lines issue.
    
    for line in lines:
        # All Lines with sequence IDs start with '>'
        if line.startswith('>'):  
            # Add sequences in temp to seqList and ID to IDList when new ID encountered    
            if not len(tempList) == 0:
                IDList.append(tempList[0]) # Append the sequence ID
                seqList.append(''.join(tempList[1:])) # Append everything but the ID
                tempList = [line,] # Empty List, Add new ID to it   
            # Add seq ID to list if it's the first sequence ID
            else:
                tempList = [line,]
        # Add sequences to the tempList if no ID is encountered
        else:
            tempList.append(line)
    # Handles sequence for last ID
    else:
        IDList.append(tempList[0]) # Append the sequence ID
        seqList.append(''.join(tempList[1:])) # Append everything but the ID

assert len(seqList) == 1


# Need to accomplish the 2 below in a loop. Should assign ID to seq and save everything in a dictionary. IDs = key
    # Value = list with coord and protein sequence 

# TODO: need to run longestNonCoding to find min threshold for id
            
# TODO: Use output from longestNonCoding to set conservative threshold for min length needed for ORF to be real 
    # and run gene finder to find real orf.

#res = geneFinder(seqList[0], 20)
#print(res)