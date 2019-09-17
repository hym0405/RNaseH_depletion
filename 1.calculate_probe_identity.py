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
def readProbeTsv(path):
    f = open(path,"rb")
    data = f.readlines()
    f.close()
    Pool = {}
    for each in data[1:]:
        each = each.rstrip()
        tmp = each.split("\t")
        Pool[tmp[0]] = []
    for each in data[1:]:
        each = each.rstrip()
        tmp = each.split("\t")
        Pool[tmp[0]].append([tmp[1], tmp[2]])
    return Pool
def convertProbeToRNA(pool):
    outPool = {}
    for eachSet in pool.keys():
        probeSet = pool[eachSet]
        tmpPool = []
        for eachProbe in probeSet:
            tmpID = eachProbe[0]
            tmpSeq = eachProbe[1]
            tmp = tmpID.split("_")
            tmpOrder = int(tmp[-1])
            tmpPool.append([tmpOrder, tmpSeq])
        tmpPool.sort()
        mergedSeq = "".join([e[1] for e in tmpPool])
        rawRNA = reverseComplement(mergedSeq)
        tmp = eachSet.split("_")
        RNAtype = tmp[-1]
        probeSet.reverse()
        outPool[RNAtype] = [eachSet, rawRNA, probeSet]
    return outPool
def reverseComplement(seq):
    basePool = {"A":"T","T":"A","C":"G","G":"C"}
    out_seq = ""
    for e in seq:
        out_seq = basePool[e] + out_seq
    return out_seq
def calculateMis(ref, sub):
    count = 0
    for i in range(len(ref)):
        if ref[i] == sub[i]:
            pass
        else:
            count += 1
    return count
def muscleAlignment(probeRNA_ID, probeRNA_seq, targetRNA_ID, targetRNA_seq):
    f = open("./tmp.fa","w")
    f.writelines([">" + probeRNA_ID + os.linesep])
    f.writelines([probeRNA_seq + os.linesep])
    f.writelines([">" + targetRNA_ID + os.linesep])
    f.writelines([targetRNA_seq + os.linesep])
    f.close()
    os.system("chmod +x " + muscle_exe_PATH)
    os.system(muscle_exe_PATH + " -in ./tmp.fa -out ./tmp.alignment.fa")
    probe_target_alignment_list = readFASTA("./tmp.alignment.fa")
    os.system("rm -f ./tmp.fa ./tmp.alignment.fa")
    return probe_target_alignment_list
def findProbeMap(alignment_seq, original_seq, probe_set):
    probeLength_list = [len(e) for e in probe_set]
    probeIndex_original = []
    for i in range(len(probeLength_list)):
        for j in range(probeLength_list[i]):
            probeIndex_original.append(i)
    probeIndex_alignment = []
    probeLength_alignment = []
    for e in probeLength_list:
        probeLength_alignment.append(0)
    index = 0
    for e in alignment_seq:
        probeIndex_alignment.append(probeIndex_original[index])
        probeLength_alignment[probeIndex_original[index]] += 1
        if e != "-":
            index += 1  
    return probeIndex_alignment, probeLength_alignment
def findTargetMap(alignment_seq, location):
    count = 0
    for i in range(location):
        if alignment_seq[i] != "-":
            count += 1
    return count

#######################
### Input Arguments ###
#######################


parser = argparse.ArgumentParser(description = "Calculate probe identity to new rRNA sequences to evaluate the ability of pools to be applied to different sequences")
parser.add_argument("-t", "--target", type = str,
                    help="Path to target rRNA sequences. All rRNA sequences should be labelled as [SampleID]_16S and [SampleID]_23S in FASTA format")
parser.add_argument("-p", "--probe", type = str, 
                    help="Path to probe sequences to be evaluated")
parser.add_argument("-o", "--output_prefix", type = str,  
                    help="Prefix of output probe identity file. Results of probe identity for different rRNA sequences will be saved as individual files labelled as [output_prefix].[rRNA_Label].tsv")
parser.add_argument("-m", "--muscle_path", type = str, default = "./bin/muscle",
					help="Path to executable file of muscle [default: ./bin/muscle]")
args = parser.parse_args()



##Path to target rRNA sequences. All rRNA sequences should be labelled as [SampleID]_16S 
##and [SampleID]_23S in FASTA format.
##Examples can be found in ./data/rRNA_sequence.dorei.fa
#input_target_rRNA_PATH = "./data/rRNA_sequence.aerofaciens.fa"
input_target_rRNA_PATH = args.target

##Path to probe sequences
##Examples can be found in ./output/rRNA_probe.dorei.tsv
#input_probe_PATH = "./output/rRNA_probe.dorei.tsv"
input_probe_PATH = args.probe

##Prefix of output probe identity file.
##Results of probe identity for different rRNA sequences will be saved as individual files
##labelled as [output_probe_identity_Prefix].[rRNA_Label].tsv
#output_probe_identity_Prefix = "./output/probeIdentity.probe_dorei"
output_probe_identity_Prefix = args.output_prefix


##The Path to executable file of muscle
#muscle_exe_PATH = "./bin/muscle"
muscle_exe_PATH = args.muscle_path

###############################################
### Program for Probe Indentity Calculation ###
###############################################

##Read probe sequences
probe_pool = readProbeTsv(input_probe_PATH)

##Merge probe sequences together and convert them back to original rRNA sequences
probeRNA_pool = convertProbeToRNA(probe_pool)

##Read target rRNA sequences
targetrRNA_list = readFASTA(input_target_rRNA_PATH)

for eachTargetRNA in targetrRNA_list:
    targetRNA_ID = eachTargetRNA[0]
    targetRNA_seq = eachTargetRNA[1]
    tmp = targetRNA_ID.split("_")
    targetRNA_type = tmp[-1]
    probeRNA_set = probeRNA_pool[targetRNA_type]
    probeRNA_ID = probeRNA_set[0]
    probeRNA_seq = probeRNA_set[1]
    probeSet_seq = probeRNA_set[2]
    
    ##Run Muscle to get the alignment of probe rRNA and target rRNA
    probe_target_alignment_list = muscleAlignment(probeRNA_ID, probeRNA_seq, targetRNA_ID, targetRNA_seq)
    probeRNA_alignedSeq = probe_target_alignment_list[0][1]
    targetRNA_alignedSeq = probe_target_alignment_list[1][1]
    length_of_alignment = len(targetRNA_alignedSeq)
    
    ##Calculate total number of mismatches
    num_of_total_mismatches = calculateMis(probeRNA_alignedSeq, targetRNA_alignedSeq)
    
    ##Calculate the number of mismatches for each individual probe 
    probe_alignedIndex, probe_alignedLength = findProbeMap(probeRNA_alignedSeq, probeRNA_seq, [e[1] for e in probeSet_seq])
    probe_mismatch_pool = {}
    for i in range(len(probeSet_seq)):
        probe_mismatch_pool[i] = 0
    for i in range(len(targetRNA_alignedSeq)):
        probeBase = probeRNA_alignedSeq[i]
        targetBase = targetRNA_alignedSeq[i]
        probeIndex = probe_alignedIndex[i]
        if probeBase != targetBase:
            probe_mismatch_pool[probeIndex] += 1
            
    ##Write the results of probe identity to [output_probe_identity_Prefix] for each target rRNA 
    f = open(output_probe_identity_Prefix + "." + targetRNA_ID + ".tsv","w")
    f.writelines(["## Target rRNA:" + targetRNA_ID + os.linesep])
    f.writelines(["## Probe set designed for: " + probeRNA_ID + os.linesep])
    f.writelines(["## Total length of target rRNA " + targetRNA_ID + ": " + str(len(targetRNA_seq)) + os.linesep])
    f.writelines(["## Total length of probe-target alignment: " + str(length_of_alignment) + os.linesep])
    f.writelines(["## Number of mismatches in probe-target alignment: " + str(num_of_total_mismatches) + os.linesep])
    f.writelines(["#target_ID\ttarget_start\ttarget_end\tprobe_ID\tlength_alignment\tnum_of_mismatches\tratio" + os.linesep])
    for i in range(len(probeSet_seq)):
        target_ID = targetRNA_ID
        probe_ID = probeSet_seq[i][0]
        length_alignment = probe_alignedLength[i]
        target_start = findTargetMap(targetRNA_alignedSeq, sum(probe_alignedLength[:i]) + 1)
        target_end = findTargetMap(targetRNA_alignedSeq, sum(probe_alignedLength[:i + 1]))
        num_of_mismatches = probe_mismatch_pool[i]
        ratio = 1.0 * num_of_mismatches / length_alignment
        toWrite = [target_ID, target_start, target_end, probe_ID, length_alignment, num_of_mismatches, ratio]
        f.writelines(["\t".join([str(e) for e in toWrite]) + os.linesep])
    f.close()
