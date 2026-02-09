#!/usr/bin/env python

import sys
import os

r1_path = sys.argv[1]
r2_path = sys.argv[2]
output_path = sys.argv[3]

def read_fasta_pairs(file_path):
    with open(file_path) as f:
        lines = f.readlines()
    return [(lines[i], lines[i+1]) for i in range(0, len(lines), 2)]

r1_pairs = read_fasta_pairs(r1_path)

r2_exists = os.path.isfile(r2_path)
if r2_exists:
    r2_pairs = read_fasta_pairs(r2_path)

output_lines = []
for i in range(len(r1_pairs)):
    output_lines.extend(r1_pairs[i])
    if r2_exists:
    	output_lines.extend(r2_pairs[i])

with open(output_path, 'w') as out:
    out.writelines(output_lines)
