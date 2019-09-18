#!/usr/bin/env python
import os 
import sys
import argparse
def tsv2Fasta(In, Out):
	f = open(In,"rb")
	data = f.readlines()
	f.close()
	f = open(Out,"w")
	for each in data:
		if each[:-1] == "rRNA_label\tprobe_ID\tprobe_sequence":
			continue
		each = each[:-1]
		tmp = each.split("\t")
		f.writelines([">" + tmp[1] + os.linesep])
		f.writelines([tmp[2] + os.linesep])
	f.close()
def processBLASToutput(inputF, reference, outputF):
	f = open(reference,"rb")
	data = f.readlines()
	f.close()
	Pool = []
	for each in data:
		if each[-2:] == "\r\n":
			Pool.append(each[:-2])
		elif each[-1] == "\n":
			Pool.append(each[:-1])
	f = open(inputF, "rb")
	data = f.readlines()
	f.close()
	f = open(outputF,"w")
	f.writelines("probeID\ttranscript(off-target)" + os.linesep)
	for each in data:
		each = each[:-1]
		if each[0] == "#":
			continue
		else:
			tmp = each.split("\t")
			if tmp[1] not in Pool:
				f.writelines([tmp[0] + "\t" + tmp[1] + os.linesep])
	f.close()
def processBURSToutput(inputF, reference, outputF):
	num_of_mismatch = 8
	f = open(reference,"rb")
	data = f.readlines()
	f.close()
	Pool = []
	for each in data:
		if each[-2:] == "\r\n":
			Pool.append(each[:-2])
		elif each[-1] == "\n":
			Pool.append(each[:-1])
	f = open(inputF, "rb")
	data = f.readlines()
	f.close()
	f = open(outputF,"w")
	f.writelines("probeID\ttranscript(off-target)" + os.linesep)
	for each in data:
		each = each[:-1]
		if each[0] == "#":
			continue
		else:
			tmp = each.split("\t")
			tmp2 = tmp[1].split(" ")
			if tmp2[0] not in Pool and int(tmp[4]) <= num_of_mismatch:
				f.writelines([tmp[0] + "\t" + tmp2[0] + os.linesep])
	f.close()

parser = argparse.ArgumentParser(description = "Predict potential off-targets for probe libraries")
parser.add_argument("-t", "--transcript", type = str,
                    help="Path to transcript sequences. All transcript sequences should be saved in FASTA format")
parser.add_argument("-r", "--rRNA", type = str,
					help="Path to list of rRNA transcript IDs")
parser.add_argument("-p", "--probe", type = str,
                    help="Path to probe sequences to be evaluated. Probe sequences can be saved in either TSV or FASTA format (should be specified in probe format)")
parser.add_argument("-pf", "--probe_format", choices=["TSV", "FASTA"], default = "TSV", 
					help="Format of probe sequences, either TSV or FASTA [default: TSV]")
parser.add_argument("-o", "--output_prefix", type = str,
                    help="Prefix of output predicted off-targets file. Results of BLASTN and BURST will be saved in [output_prefix].BLAST.tsv and [output_prefix].BURST.tsv")
parser.add_argument("-mb", "--makeblastdb_path", type = str, default = "./bin/makeblastdb",
                    help="Path to executable file of makeblastdb (NCBI-BLAST) [default: ./bin/makeblastdb]")
parser.add_argument("-bn", "--blastn_path", type = str, default = "./bin/blastn",
                    help="Path to executable file of blastn (NCBI-BLAST) [default: ./bin/blastn]")
parser.add_argument("-br", "--burst_path", type = str, default = "./bin/burst",
                    help="Path to executable file of burst [default: ./bin/burst]")
args = parser.parse_args()	

path_transcript = args.transcript
path_rRNAref = args.rRNA
path_probe = args.probe
path_outputP = args.output_prefix
format_probe = args.probe_format
bin_makeblastdb = args.makeblastdb_path
bin_blastn = args.blastn_path
bin_burst = args.burst_path

os.system("chmod +x " + bin_makeblastdb)
os.system("chmod +x " + bin_blastn)
os.system("chmod +x " + bin_burst)
os.system("mkdir -p ./.tmp_offtarget_predict_RNaseH")

if format_probe == "FASTA":
	os.system("cp " + path_probe + " ./.tmp_offtarget_predict_RNaseH/probe.fa")
else:
	tsv2Fasta(path_probe, "./.tmp_offtarget_predict_RNaseH/probe.fa")

os.system("cp " + path_transcript + " ./.tmp_offtarget_predict_RNaseH/target.ffn")
os.system(bin_makeblastdb + " -in ./.tmp_offtarget_predict_RNaseH/target.ffn -dbtype nucl -parse_seqids " + \
		"-out ./.tmp_offtarget_predict_RNaseH/target")
os.system(bin_blastn + " -query ./.tmp_offtarget_predict_RNaseH/probe.fa -db ./.tmp_offtarget_predict_RNaseH/target " + \
		"-outfmt \"7 qaccver saccver pident evalue bitscore qcovs length mismatch gapopen\" " + \
		"-num_threads 4 -out ./.tmp_offtarget_predict_RNaseH/target.blast.out")
processBLASToutput("./.tmp_offtarget_predict_RNaseH/target.blast.out", \
					path_rRNAref, path_outputP + ".BLAST.tsv")

os.system(bin_burst + " --references ./.tmp_offtarget_predict_RNaseH/target.ffn --queries " + \
		"./.tmp_offtarget_predict_RNaseH/probe.fa -fr -i 0.80 -m FORAGE --threads 4 --output " + \
		"./.tmp_offtarget_predict_RNaseH/target.burst.out")
processBURSToutput("./.tmp_offtarget_predict_RNaseH/target.burst.out", \
					 path_rRNAref, path_outputP + ".BURST.tsv")
os.system("rm -rf ./.tmp_offtarget_predict_RNaseH")
