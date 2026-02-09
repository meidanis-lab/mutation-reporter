#!/usr/bin/env python


''' Program to find mutations

    Developed at the Boldrini Center, Brazil, in 2023.
'''

############################################################
### imports

import os
import glob
import sys
import argparse
import collections
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Blast import NCBIXML

############################################################
### constants

MSGS = False

############################################################
### classes

class Mutation:

    def __init__(self, transcript_id, gene, sbjct_pos, pos, query_pos, count):
        self.transcript_ids = [transcript_id]
        self.gene = gene
        self.sbjct_pos = sbjct_pos
        self.pos = pos
        self.query_pos = query_pos
        self.count = count


############################################################
### routines

############################################################
### routine depth_index

def depth_index(vector, first, last, item):
    '''
    routine depth_index

    parameters:
    vector: integer vector sorted in strictly increasing order
    first, last: indexes for vector, such that first <= last+1
    item: search key such that vector[first] <= item

    return value:
    largest index i such that vector[i] <= item
    '''

    if last >= first:
        middle = (first + last) // 2
        if vector[middle] <= item:
            return depth_index(vector, middle + 1, last, item)
        else:
            return depth_index(vector, first, middle - 1, item)
    else:
        return last


############################################################
### routine update_mutations

def update_mutations(mutations, conserv, transcript_id, hsp, gene):
    '''
    routine update_mutations

    parameters:
    mutations: mutation dictionary to be updated
    conserv: caracter representing mutation type:
             ' ': amino acid mismatch
             '+': amino acid conserved substitution
    transcript_id: id of transcript as reported by Illumina fastq
    hsp: high scoring pair BLAST structure
    gene: gene where mutation was found

    return value:
    none, mutation dictionary updated
    '''

    pos = -1
    while True:
        pos = hsp.match.find(conserv, pos + 1)
        if pos < 0:
            break
        protmut = hsp.sbjct[pos] + str(pos + hsp.sbjct_start) + hsp.query[pos]
        if protmut in mutations:
            mutations[protmut].count += 1
            mutations[protmut].transcript_ids.append(transcript_id)
        else:
            mutations[protmut] = Mutation(
                transcript_id,
                gene,
                hsp.sbjct[pos],
                str(pos + hsp.sbjct_start),
                hsp.query[pos],
                1)

############################################################
### routine construct_transcripts

def construct_transcripts(parse, msgs, min_align_len, min_perc_ident):
    '''
    routine construct_transcripts

    parameters:
    parse: parsing object for BLAST xml output
    msgs: debug flag
    min_align_len: minimum length for an alignment to be included
    min_perc_ident: minimum percent identity for an alignment to be included

    return value:
    list_start: list of starting positions in the reference of all included alignments
    list_end: list of ending positions in the reference of all included alignments
    mutations: dictionary of mutations found (see class Mutation)
    transcript_list: list of transcripts; each transcript has two parts; each part has a start and end
    count_hits: total number of alignments in parse
    count_good: number of included alignments

    note: list_start and list_end usually contain many repetitions
    note: usually, each record has just one alignment and each alignment has just one hsp
    '''

    count_hits = 0
    ## counting good hits (count_good)
    ## a good hit is one that satisfies the criteria
    ## for length and percentual identity
    count_good = 0
    list_start = []
    list_end = []
    list_transcript_id = []
    transcript_list = []
    waiting = False
    mutations = {}
    for record in parse:
        count_hits += 1
        if record.alignments:
            count_align = 0
            for align in record.alignments:
                count_align += 1
                if count_align > 1:
                    if msgs:
                        print('More than one hit for ' + record.query.split(' ')[1])
                    break
                count_hsp = 0
                for hsp in align.hsps:
                    count_hsp += 1
                    if count_hsp > 1:
                        if msgs:
                            print('More than one hsp for ' + record.query)
                    percent_ident = 100 * hsp.identities / hsp.align_length

                    if msgs:
                       print(record.query_id)
                       print(record.query)
                       print(hsp.query)
                       print(hsp.match)
                       print(hsp.sbjct)
                       print(hsp.expect)
                       print(percent_ident)

                    if hsp.align_length >= int(min_align_len):
                        if percent_ident > int(min_perc_ident):

                            count_good += 1

                            list_start.append(hsp.sbjct_start)
                            list_end.append(hsp.sbjct_end + 1)

                            transcript_id = (':'.join(record.query.split()[0].split(':')[-2:]))
                            list_transcript_id.append(transcript_id)

                            if len(list_transcript_id) >= 2 and list_transcript_id[-1] == list_transcript_id[-2]:
                                transcript = (list_start[-2],
                                          list_end[-2],
                                          hsp.sbjct_start,
                                          hsp.sbjct_end)
                                transcript_list.append(transcript)
                                waiting = False
                            else:
                                if waiting:
                                    transcript = (list_start[-2],
                                          list_end[-2],
                                          list_start[-2],
                                          list_end[-2])
                                    transcript_list.append(transcript)
                                    waiting = False
                                else:
                                    waiting = True

                            update_mutations(mutations, ' ', transcript_id, hsp, align.hit_id)
                            update_mutations(mutations, '+', transcript_id, hsp, align.hit_id)

    return list_start, list_end, mutations, transcript_list, count_hits, count_good

############################################################
### routine cummulative_counts

def cummulative_counts(list_pos):
    '''
    routine cummulative_counts

    parameters:
    list_pos: a list of positions, usually with lots or repetitions

    return value:
    pos: vector of sorted positions found in list_pos
    cumm: number of items in list that are <= the corresponding pos value
    '''

    counter = collections.Counter(list_pos)
    dictionary = dict(sorted(counter.items()))

    acc = 0
    pos = [0]
    cumm = [0]
    for key, val in dictionary.items():
        acc += val
        pos.append(key)
        cumm.append(acc)

    return pos, cumm

############################################################
### routine print_individual_mutations

def print_individual_mutations(mutations, pos_start, pos_end,
                               value_start, value_end,
                               min_depth, min_vaf):
    '''
    routine print_individual_mutations

    prints all mutations found that comply with the VAF and depth
    restrictions
    
    parameters:
    mutations: dictionary of mutations
    pos_start: sorted vector of transcript start positions
    pos_end: sorted vector of transcript end positions
    value_start: number of transcripts for the corresponding start position
    value_end: number of transcripts for the corresponding end position
    min_depth: minimum depth needed to report a mutation
    min_vaf: minimum VAF needed to report a mutation
    
    return value:
    valid_mutations: list of mutations that satisfy the criteria

    collateral effect:
    prints a table with mutations that satisfy the criteria
    '''

    valid_mutations = []

    print('Gene', '\t', 'Prot Mut', '\t', 'Mut Reads', '\t', 'WT Reads', '\t', 'VAF')

    for m in mutations:
        start_depth = value_start[(depth_index(pos_start, 0, len(pos_start) - 1, int(m[1:-1])))]
        end_depth = value_end[(depth_index(pos_end, 0, len(pos_end) - 1, int(m[1:-1])))]
        depth = start_depth - end_depth
        if depth > int(min_depth):
            VAF = 100 * mutations[m].count / depth
            if VAF > int(min_vaf):
                valid_mutations.append((m, mutations[m]))
                print(mutations[m].gene, '\t', m, '\t', mutations[m].count,'\t', depth - mutations[m].count, '\t', round(VAF, 2), '%')

    return valid_mutations

############################################################
### routine print_mutation_pairs

def print_mutation_pairs(valid_mutations, transcript_list, min_depth):
    '''
    routine print_mutation_pairs

    prints all pairs of same-transcript mutations
    with depth larger than min_depth
    
    parameters:
    valid_mutations: list of reportable mutations
    transcript_list: list of transcripts, each given by
        the coordinates of their two reads
    min_depth: threshold for reporting a mutation pair
    
    return value:
    none

    collateral effect:
    prints a table with mutation pairs that satisfy the criteria
    '''

    print()
    print('Gene', '\t', 'Prot Mut 1', '\t', 'Prot Mut 2', '\t', 'Double Mutated Transcripts', '\t', 'Total Transcripts', '\t', 'VAF')
    nval = len(valid_mutations)
    reads_valid_mutations_i = []
    reads_valid_mutations_j = []

    for i in range(nval):
        posi_i = int(valid_mutations[i][0][1:-1])
        for j in range(i+1, nval):
            depth = 0
            posi_j = int(valid_mutations[j][0][1:-1])
            reads_valid_mutations_i = set(valid_mutations[i][1].transcript_ids)
            reads_valid_mutations_j = set(valid_mutations[j][1].transcript_ids)
            numerator = len(set(reads_valid_mutations_i) & set(reads_valid_mutations_j))
            for transcript in transcript_list:
                a1, a2, b1, b2 = transcript
                if ((min(a1, a2) <= posi_i <= max(a1,a2)) or (min(b1, b2) <= posi_i <= max(b1,b2))) and ((min(a1, a2) <= posi_j <= max(a1,a2)) or (min(b1, b2) <= posi_j <= max(b1,b2))):
                    depth += 1
            if depth > int(min_depth):
                VAF = 100 * numerator / depth
                print(valid_mutations[i][1].gene, '\t', valid_mutations[i][0], '\t', valid_mutations[j][0], '\t', numerator, '\t', depth, '\t', round(VAF, 2), '%')


############################################################
### routine that parses blast XML output
###
### xml: file generated by BLAST

def parser(min_align_len, min_perc_ident, min_depth, min_vaf, xml, msgs):
    parse = NCBIXML.parse(open(xml))
    mut = []

    list_start, list_end, mutations, transcript_list, count_hits, count_good = construct_transcripts(parse, msgs, min_align_len, min_perc_ident)

    ## Compute cummulative counts for start and end positions
    pos_start, value_start = cummulative_counts(list_start)
    pos_end,   value_end   = cummulative_counts(list_end)

    valid_mutations = print_individual_mutations(
        mutations, pos_start, pos_end, value_start, value_end,
        min_depth, min_vaf)
    print_mutation_pairs(valid_mutations, transcript_list, min_depth)

    return(count_hits, count_good)

def main(MSGS):
    MIN_ALIGN_LEN = sys.argv[1]
    MIN_PERC_IDENT = sys.argv[2]
    MIN_DEPTH = sys.argv[3]
    MIN_VAF = sys.argv[4]
    XML = sys.argv[5]
    if len(sys.argv) >= 7:
        MSGS = True

    c_h, c_g = parser(MIN_ALIGN_LEN, MIN_PERC_IDENT, MIN_DEPTH,
                      MIN_VAF, XML, MSGS)

if __name__ == "__main__":
    main(MSGS)

