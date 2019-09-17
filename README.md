# Scalable and cost-effective ribonuclease-based rRNA depletion for bacterial transcriptomics

Tools for probe design and evaluation:

* **0.design_probe.py**: design probe libraries for target 16S and 23S rRNA sequence

* **1.calculate_probe_identity.py**: calculate probe identity to various different 16S and 23S sequences to evaluate the ability of pools to be applied to different sequences

* **2.predict_probe_offtarget.sh**: predict potential off-targets for probe libraries

## Dependencies

* Python 2.7, Jupyter 4.3.0 
	- panda
	- numpy
	- argparse
	- _NB: Above libraries are bundled together in the [Anaconda distribution](https://www.continuum.io/downloads)_


* [Muscle: MUltiple Sequence Comparison by Log-Expectation](https://www.drive5.com/muscle/)
	- **Required for probe identity calculation only**
	- Executable file of Muscle that compatible with your operating system should be put into ./bin or other place specified in 1.calculate_probe_identity.ipynb
	
* [NCBI BLAST+ Executables 2.9.0](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
	- **Required for probe off-targets prediction only**
	- Executable file of makeblastdb and blastn that compatible with your operating system should be put into ./bin or other place specified in 2.predict_probe_offtarget.sh
	
* [BURST v0.99.8 DB15](https://github.com/knights-lab/BURST)
	- **Required for probe off-targets prediction only**
	- Executable file of burst that compatible with your operating system should be put into ./bin or other place specified in 2.predict_probe_offtarget.sh


## Design probe libraries for 16S and 23S rRNA sequence

### Description
```
usage: 0.design_probe.py [-h] [-i INPUT] [-o OUTPUT] [-l LENGTH]

Design probe libraries for bacterial rRNA depletion

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to rRNA sequences. All rRNA sequences should be
                        labelled as [SampleID]_16S and [SampleID]_23S in FASTA
                        format
  -o OUTPUT, --output OUTPUT
                        Path to output probe file. Probe sequences will be
                        saved as a tab-delimited table
  -l LENGTH, --length LENGTH
                        Length of probes [default: 50]
```
### Input format
****[Important] Avoid underline in sample IDs****

**rRNA sequence:** 16S and 23S rRNA sequence in FASTA format and all rRNA sequences should be labelled as [SampleID]_16S and [SampleID]_23S

**** We used [Prokka](https://github.com/tseemann/prokka) to predict 16S and 23S rRNA sequences of bacterial species****

****[example: ./data/rRNA_sequence/rRNA_sequence.dorei.fa]****

```
>dorei_16S
AGAGTTTGATCCTGGCTC...
...
>dorei_23S
GAAAGTAAAGAAGGGCGC...
...
```
### Output format
**** Probe sequences will be saved as a tab-delimited table****

****[example: ./output/rRNA_probe.dorei.tsv]****
```
rRNA_label      probe_ID        probe_sequence
dorei_16S       dorei_16S_0     AGGTGTTCCAGCCGC...
dorei_16S       dorei_16S_1     GTTTTACCCTAGGGC...
dorei_16S       dorei_16S_2     TCCCATGGCTTGACG...
...		...		...
dorei_23S       dorei_23S_0     TAAGGAAAGTGGACG...
dorei_23S       dorei_23S_1     CAACGTCGTAGTCTA...
dorei_23S       dorei_23S_2     TCGTACTTAGATGCT...
...
```

### Example
```
chmod +x ./0.design_probe.py
./0.design_probe.py -i ./data/rRNA_sequence/rRNA_sequence.dorei.fa \
		-o ./output/rRNA_probe.dorei.tsv \
		-l 50
```





