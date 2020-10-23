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

