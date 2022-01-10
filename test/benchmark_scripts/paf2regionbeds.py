#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 15:12:09 2022

@author: Jordan
"""

import sys
import os

if __name__ == "__main__":
    
    paf_name = sys.argv[1]
    slop = int(sys.argv[2])
    min_mapq = int(sys.argv[3])
    max_len = int(sys.argv[4])
    query_bed_name = sys.argv[5]
    ref_bed_name = sys.argv[6]
    
    assert not os.path.exists(query_bed_name)
    assert not os.path.exists(ref_bed_name)
    
    contig_strand_count = {}
    for line in open(paf_name):
        tokens = line.strip().split()
        qname = tokens[0]
        rev = tokens[4] == "-"
        rname = tokens[5]
        alnlen = int(tokens[3]) - int(tokens[2])
        
        if qname not in contig_strand_count:
            contig_strand_count[qname] = {}
        key = (rname, rev)
        if key not in contig_strand_count[qname]:
            contig_strand_count[qname][key] = alnlen
        else:
            contig_strand_count[qname][key] += alnlen
            
    
    contig_to_strand = {}
    for qname in contig_strand_count:
        contig_to_strand[qname] = max(contig_strand_count[qname], key = lambda k: contig_strand_count[qname][k])
    
    for qname in contig_to_strand:
        rname, rev = contig_to_strand[qname]
        if rev:
            strand = "-"
        else:
            strand = "+"
        print("{} attributed to {}{}".format(qname, rname, strand), file = sys.stderr)
    
    curr_qname = None
    curr_rbegin = None
    curr_rend = None
    curr_qbegin = None
    curr_qend = None
    
    query_bed = open(query_bed_name, "w")
    ref_bed = open(ref_bed_name, "w")
    
    def emit(curr_qname, curr_qbegin, curr_qend, curr_rbegin, curr_rend):
        if curr_qname is not None:
            print("{}\t{}\t{}\t.\t.\t+".format(curr_qname, curr_qbegin, curr_qend), file = query_bed)
            curr_rname, rev = contig_to_strand[qname]
            if rev:
                strand = "-"
            else:
                strand = "+"
            print("{}\t{}\t{}\t.\t.\t{}".format(curr_rname, curr_rbegin, curr_rend, strand), file = ref_bed)
        return qname, qbegin, qend, rbegin, rend
    
    for line in open(paf_name):
        tokens = line.strip().split()
        qname = tokens[0]
        qbegin = int(tokens[2])
        qend = int(tokens[3])
        rev = tokens[4] == "-"
        rname = tokens[5]
        rbegin = int(tokens[7])
        rend = int(tokens[8])
        mapq = int(tokens[11])
        
        # low mapq or mismatched mapping, skip entirely
        if mapq < min_mapq or (rname, rev) != contig_to_strand[qname]:
            continue
            
        # first block
        if curr_rbegin is None:
            curr_qname, curr_qbegin, curr_qend, curr_rbegin, curr_rend = qname, qbegin, qend, rbegin, rend
            continue
        
        qdiff = abs(qbegin - curr_qend - 1)
        if rev:
            rdiff = abs(curr_rbegin - rend - 1)
        else:
            rdiff = abs(rbegin - curr_rend - 1)
        
        next_qlen = qend - curr_qbegin + 1
        next_rlen = max(curr_rend, rend) - min(curr_rbegin, rbegin) + 1
        
        # if too long, too separated, or between contigs, emit the current and reset
        if qdiff > slop or rdiff > slop or qname != curr_qname or next_qlen > max_len or next_rlen > max_len:
            emit(curr_qname, curr_qbegin, curr_qend, curr_rbegin, curr_rend)
            curr_qname, curr_qbegin, curr_qend, curr_rbegin, curr_rend = qname, qbegin, qend, rbegin, rend
            continue
        
        curr_qname = qname
        curr_qbegin = min(qbegin, curr_qbegin)
        curr_qend = max(qend, curr_qend)
        curr_rbegin = min(rbegin, curr_rbegin)
        curr_rend = max(rend, curr_rend)
        
    # emit the final
    emit(curr_qname, curr_qbegin, curr_qend, curr_rbegin, curr_rend)
        
            
            
        
        
        
        
    
    