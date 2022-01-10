#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 21:42:47 2022

@author: Jordan
"""

import sys
import os

def parse_cigar(cigar):
    
    qlen = 0
    rlen = 0
    
    sclip_start = 0
    sclip_end = 0
    j = 0
    for i in range(len(cigar)):
        if not cigar[i].isdigit():
            olen = int(cigar[j:i])
            op = cigar[i]
            if op == "S" or op == "H":
                if i == len(cigar) - 1:
                    sclip_end += olen
                else:
                    sclip_start += olen
            elif op == "M"  or op == "=":
                qlen += olen
                rlen += olen
            elif op == "X":
                qlen += olen
                rlen += olen
            elif op == "I":
                qlen += olen
            elif op == "D" or op == "N":
                rlen += olen
            j = i + 1
        
    return qlen, rlen, sclip_start, sclip_end

if __name__ == "__main__":
    
    min_length = int(sys.argv[1])
    min_mapq = int(sys.argv[2])
    fasta_out_dir = sys.argv[3]
    summary_file = sys.argv[4]
    
    assert(not os.path.exists(fasta_out_dir))
    assert(not os.path.exists(summary_file))
    
    os.mkdir(fasta_out_dir)
    
    has_dupes = set()
    table = {}
        
    i = 0
    for line in sys.stdin:
        
        i += 1
        if i % 10000 == 0:
            print("processed " + str(i) + " reads", file = sys.stderr)
        tokens = line.strip().split()
        seq = tokens[9]
        if len(seq) < min_length:
            continue
        name = tokens[0]
        fa_name = os.path.join(fasta_out_dir, name + ".fa")
        if os.path.exists(fa_name):
            has_dupes.add(name)
            del table[name]
            os.remove(fa_name)
            continue
        if name in has_dupes:
            continue
        
        mapq = int(tokens[4])
        if mapq < min_mapq:
            has_dupes.add(name)
            continue
        
        qaln, raln, scs, sce = parse_cigar(tokens[5])
        
        flags = int(tokens[1])
        rev = bool(flags & 16)
        rname = tokens[2]
        
        pos_begin = int(tokens[3])
        pos_end = pos_begin + raln - 1
            
        strand = "+"
        if rev:
            strand = "-"
        
        table[name] = "\t".join(map(str, [rname, strand, pos_begin, pos_end]))
        
        with open(fa_name, "w") as f:
            print(">" + name + " " + str(scs) + ":" + str(len(seq) - sce) + " of " + str(len(seq)), file = f)
            print(seq[scs:len(seq)-sce], file = f)
        
        
    with open(summary_file, "w") as f:
        for name in sorted(table):
            print(name + "\t" + table[name], file = f)
        

