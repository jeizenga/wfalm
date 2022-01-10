#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 16:40:09 2022

@author: Jordan
"""

import sys
import os
import subprocess


def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

def rev_comp_fasta(fa):
    seq = ""
    header = ""
    for line in open(fa):
        if line.startswith(">"):
            assert(header == "")
            header = line.strip()
        else:
            seq += line.strip()
            
    with open(fa, "w") as f:
        print(header, file = f)
        print(reverse_complement(seq), file = f)

if __name__ == "__main__":
    
    bed_name = sys.argv[1]
    fasta_name = sys.argv[2]
    prefix = sys.argv[3]
    
    i = 1
    
    for line in open(bed_name):
        tokens = line.strip().split()
        seq = tokens[0]
        b = int(tokens[1]) + 1
        e = int(tokens[2])
        rev = tokens[5] == "-"
        
        fa_out = "{}{}.fa".format(prefix, i)
        cmd = "samtools faidx {} {}:{}-{} > {}".format(fasta_name, seq, b, e, fa_out) 
        subprocess.check_call(cmd, shell = True)
        
        if rev:
            rev_comp_fasta(fa_out)
        
        i += 1
    