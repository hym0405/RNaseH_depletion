# Scalable and cost-effective ribonuclease-based rRNA depletion for bacterial transcriptomics

Jupyter Notebook for probe design and evaluation:

* **0.design_probe.ipynb**: design probe libraries for target 16S and 23S rRNA sequence

* **1.calculate_probe_identity.ipynb**: calculate probe identity to various different 16S and 23S sequences to evaluate the ability of pools to be applied to different sequences

* **2.predict_probe_offtarget.sh**: predict potential off-targets for probe libraries

## Dependencies

* Python 2.7, Jupyter 4.3.0 
	- panda
	- numpy
	- _NB: Above libraries are bundled together in the [Anaconda distribution](https://www.continuum.io/downloads)_

* [Muscle: MUltiple Sequence Comparison by Log-Expectation](https://www.drive5.com/muscle/)
	- **Required for "1.calculate_probe_identity.ipynb" only**
	- Executable file of Muscle that compatible with your operating system should be put into ./bin or other place specified in 1.calculate_probe_identity.ipynb
	
* [NCBI BLAST+ executables 2.9.0](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
	- **Required for "2.predict_probe_offtarget.sh" only**
	- Executable file of makeblastdb and blastn that compatible with your operating system should be put into ./bin or other place specified in "2.predict_probe_offtarget.sh"
	
* [BURST v0.99.8 DB15](https://github.com/knights-lab/BURST)
	- **Required for "2.predict_probe_offtarget.sh" only**
	- Executable file of burst that compatible with your operating system should be put into ./bin or other place specified in "2.predict_probe_offtarget.sh"
