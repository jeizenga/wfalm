#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 22:51:11 2022

@author: Jordan
"""

import os
import sys
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
    
    ref_fasta_out_dir = sys.argv[1]
    read_summary = sys.argv[2]
    mat_fasta = sys.argv[3]
    pat_fasta = sys.argv[4]
    
    assert(not os.path.exists(ref_fasta_out_dir))
    os.mkdir(ref_fasta_out_dir)
    
    ctg_to_fasta = {}
    for line in open(mat_fasta + ".fai"):
        tokens = line.strip().split()
        ctg_to_fasta[tokens[0]] = mat_fasta
    for line in open(pat_fasta + ".fai"):
        tokens = line.strip().split()
        ctg_to_fasta[tokens[0]] = pat_fasta
        
    i = 0
    for line in open(read_summary):
        i += 1
        if i % 50 == 0:
            print("extracted {} sequences".format(i))
        tokens = line.strip().split()
        qname = tokens[0]
        rname = tokens[1]
        strand = tokens[2]
        pbegin = tokens[3]
        pend = tokens[4]
        ref_fasta = ctg_to_fasta[rname]
        out = os.path.join(ref_fasta_out_dir, qname + ".fa")
        cmd = "samtools faidx {} {}:{}-{} > {}".format(ref_fasta, rname, pbegin, pend, out)
        subprocess.check_call(cmd, shell = True)
    