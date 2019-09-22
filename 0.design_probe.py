#!/usr/bin/env python
#############################
### Modules and Functions ###
#############################
import os
import sys
import argparse
import numpy as np
import pandas as pd
def readFASTA(path):
    f = open(path,"rb")
    data = f.readlines()
    f.close
    pool = []
    tmpLabel = ""
    tmpSeq = ""
    flag = 0
    for each in data:
        each = each.rstrip()
        if each[0] == ">":
            if flag == 1:
                pool.append([tmpLabel, tmpSeq])
            tmpLabel = each[1:]
            tmpSeq = ""
        else:
            tmpSeq += each
        flag = 1
    pool.append([tmpLabel, tmpSeq])
    return pool
def reverseComplement(seq):
    basePool = {"A":"T","T":"A","C":"G","G":"C"}
    out_seq = ""
    for e in seq:
        out_seq = basePool[e] + out_seq
    return out_seq
def designProbe(seq, probeLength):
    probes = []
    offset = (len(seq) % probeLength) / 2
    probes.append(seq[:offset])
    for i in np.arange(offset, len(seq), probeLength):
        probes.append(seq[i: i + probeLength])
    probes_final = []
    probes_final.append(probes[0] + probes[1])
    for i in range(len(probes) - 4):
        probes_final.append(probes[i + 2])
    if len(seq) % probeLength != 0:
        probes_final.append(probes[-2] + probes[-1])
    else:
        probes_final.append(probes[-2])
        probes_final.append(probes[-1])
    return probes_final

#######################
### Input Arguments ###
#######################

parser = argparse.ArgumentParser(description = "Design probe libraries for bacterial rRNA depletion")
parser.add_argument("-i", "--input", type = str,
                    help="Path to rRNA sequences. All rRNA sequences should be labelled as [SampleID]_16S and [SampleID]_23S in FASTA format")
parser.add_argument("-o", "--output", type = str, 
					help="Path to output probe file. Probe sequences will be saved as a tab-delimited table")
parser.add_argument("-l", "--length", type = int, default = 50,
					help="Length of probes [default: 50]")
args = parser.parse_args()

##Path to rRNA sequences. All rRNA sequences should be labelled as [SampleID]_16S 
##and [SampleID]_23S in FASTA format.
##Examples can be found in ./data/rRNA_sequence.dorei.fa

#input_rRNA_PATH = "./data/rRNA_sequence/rRNA_sequence.dorei.fa"
input_rRNA_PATH = args.input

##Path to output probe file. Probe sequences will be saved as a tab-delimited table
#output_probe_PATH = "./output/rRNA_probe.dorei.tsv"
output_probe_PATH = args.output

##Length of probes
#probe_length = 50
probe_length = args.length

################################
### Program for Probe Design ###
################################

##Read rRNA sequences
rRNA_pool = readFASTA(input_rRNA_PATH)

##Get reverse complementary sequence for each rRNA
rRNA_pool_rc = [[e[0], reverseComplement(e[1])] for e in rRNA_pool]

##Split the rRNA sequence into short oligos covering the entire length of its reverse complement
rRNA_pool_probes = [[e[0], designProbe(e[1], probe_length)] for e in rRNA_pool_rc]

##Write the output probe sequences to [output_probe_PATH]
f = open(output_probe_PATH, "w")
f.writelines(["rRNA_label\tprobe_ID\tprobe_sequence" + os.linesep])
for eachRNA in rRNA_pool_probes:
    rRNA_label = eachRNA[0]
    probe_pool = eachRNA[1]
    for i in range(len(probe_pool)):
        probe_ID = rRNA_label + "_" + str(i)
        f.writelines([rRNA_label + "\t" + probe_ID + "\t" + probe_pool[i] + os.linesep])
f.close()



