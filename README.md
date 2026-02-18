# RiboSeq-Code

Code used in TZD profiling (doi: ); data availability: (GEO: GSE319740)

1. Preprocessing - via command line - RiboSeq_preprocessing.ipynb --> RiboSeq_preprocessing.py
     Performs the trimming and alignment features before in depth data analysis
     Modify RiboSeq_preprocessing.ipynb and convert to py or modify and run RiboSeq_preprocessing.py directly
     dependencies: RiboSeq_preprocessing_fxns.py, runFastQC.sh, RiboSeq_runCutadapt.sh, runBowtieBuild.sh, RiboSeq_runBowtieAlign.sh
     inputs: demultiplexed fastq.gz files, adapter sequences (.txt), reference genome and filtering .fasta/.fa files

2. Data Analysis files are all run as jupyter notebooks
   dependencies: RiboSeq_Analysis_fxns.py, RiboSeq_Analysis_plots.py

   a. Metagene analysis - RiboSeq_MetageneAnalysis.ipynb
   b. Gene-level analysis - RiboSeq_GeneAnalysis.ipynb
   c. Codon-level analysis - RiboSeq_CodonAnalysis.ipynb


Environment requirements:
  Biopython
  Bowtie
  Cutadapt
  ...
