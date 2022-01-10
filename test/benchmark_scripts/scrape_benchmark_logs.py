#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 19:43:59 2022

@author: Jordan
"""

import sys
import re

if __name__ == "__main__":
    
    logs = sys.argv[1]
    
    log_contents = open(logs).read()
    if type(log_contents) == bytes:
        log_contents = log_contents.encode("utf-8")
    
    delim = "## BEGINNING NEW BENCHMARK ALIGNMENT ##"
    boundaries = [m.start() for m in re.finditer(delim, log_contents)]
    
    assert(len(delim) > 0)
    assert(boundaries[0] == 0)
    
    boundaries.append(len(log_contents))
    
    log_contents_by_run = [log_contents[boundaries[i]:boundaries[i+1]] for i in range(len(boundaries) - 1)]
    
    mem_line_start = "MEMORY: "
    match_line_start = "MATCH: "
    aln_line_start = "ALIGNMENT: "
    base_line_start = "BASELINE MEM: "
    aln_time_line_start = "ALIGN TIME: "
    total_mem_line_start = "Maximum resident set size"
    seq_lengths_line = "[progress] aligning sequences"
    length_regex = "'\S*' \(length (\d+)\)"
    
    # keys of low_mem, suffix_tree
    config_to_block = {(False,False):0, (False,True):1, (True, False):2, (True, True):3}
    
    data_columns = ["sm-c-baseline-mem", "sm-c-total-mem", "sm-c-time",
                    "sm-st-baseline-mem", "sm-st-total-mem", "sm-st-time",
                    "lm-c-baseline-mem", "lm-c-total-mem", "lm-c-time",
                    "lm-st-baseline-mem", "lm-st-total-mem", "lm-st-time"]
    
    # columns [baseline mem, total mem, align time] 
    
    table_entries = {}
    seq_lengths = {}
    
    for run_logs_joined in log_contents_by_run:
        
        run_logs = run_logs_joined.strip().split("\n")
        assert(len(run_logs) > 6)
        
        fasta1 = run_logs[1]
        fasta2 = run_logs[2]
        
        assert(run_logs[3].startswith(mem_line_start)) 
        assert(run_logs[4].startswith(match_line_start)) 
        assert(run_logs[5].startswith(aln_line_start))
        
        low_mem = (run_logs[3][len(mem_line_start):] == "low")
        suffix_tree = (run_logs[4][len(match_line_start):] == "suffix tree")
        #local = (run_logs[5][len(aln_line_start):] == "local")
        
        i = 6
        while not run_logs[i].startswith(base_line_start):
            i += 1
        
        baseline_mem = int(run_logs[i][len(base_line_start):].split()[0]) * 1024
        i += 1
        
        assert(run_logs[i].startswith(seq_lengths_line))
        matches = list(re.finditer(length_regex, run_logs[i]))
        assert(len(matches) == 2)
        len1 = int(matches[0].group(1))
        len2 = int(matches[1].group(1))
        seq_lengths[fasta1] = len1
        seq_lengths[fasta2] = len2
        
        i += 1
        
        while (not run_logs[i].startswith(aln_time_line_start) and 
               not run_logs[i].strip().startswith(total_mem_line_start)):
            i += 1
            
        didnt_crash = run_logs[i].startswith(aln_time_line_start)
        if didnt_crash:
            
            aln_time = int(run_logs[i][len(aln_time_line_start):].split()[0])
            
            while not run_logs[i].strip().startswith(total_mem_line_start):
                i += 1
            
            total_mem = int(run_logs[i].strip().split()[-1]) * 1024
            
        else:
            aln_time = -1
            total_mem = -1
            
        table_key = (fasta1, fasta2)
        if table_key not in table_entries:
            table_entries[table_key] = [float("nan") for i in range(len(data_columns))]
        
        row = table_entries[table_key]
        
        j = 3 * config_to_block[(low_mem, suffix_tree)]
        
        row[j] = baseline_mem
        row[j + 1] = total_mem
        row[j + 2] = aln_time
    
    srt_col = 7 # memory usage of low memory alg, direct match, probably most reliable indicator of size
    key_order = sorted(table_entries, key = lambda k: table_entries[k][srt_col])
    
    print("\t".join(["seq1", "seq2", "len1", "len2"] + data_columns) )
    
    for key in key_order:
        seq1, seq2 = key
        print("\t".join([seq1, seq2, str(seq_lengths[seq1]), str(seq_lengths[seq2])] + [str(v) for v in table_entries[key]]))
        
        
    
    