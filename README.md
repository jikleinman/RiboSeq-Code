# RiboSeq-Code

Code used in TZD profiling (doi: 10.64898/2026.02.16.706219); data availability: (GEO: GSE319740)


## 1. Environment requirements:

      download and install riboseq_env.yml as a conda environment before running code, or check version numbers from that file for incompatibilities
      
      python v3.10.12, biopython, bowtie, cutadapt, dnaio, fastQC, jupyter, matplotlib, pandas, samtools, scipy, seaborn, numpy, pysam
   


## 2. Preprocessing is run via command line 
RiboSeq_preprocessing.ipynb --> RiboSeq_preprocessing.py
   
     Performs the trimming and alignment features before in depth data analysis
   
     Modify RiboSeq_preprocessing.ipynb and convert to py or modify and run RiboSeq_preprocessing.py directly
   
     dependencies: RiboSeq_preprocessing_fxns.py, runFastQC.sh, RiboSeq_runCutadapt.sh, runBowtieBuild.sh, RiboSeq_runBowtieAlign.sh
   
     inputs: demultiplexed fastq.gz files, adapter sequences (.txt), reference genome and filtering .fasta/.fa files

   

## 3. Data Analysis files are all run as jupyter notebooks
   
      dependencies: RiboSeq_Analysis_fxns.py, RiboSeq_Analysis_plots.py
   
      a. Metagene analysis - RiboSeq_MetageneAnalysis.ipynb
      
      b. Gene-level analysis - RiboSeq_GeneAnalysis.ipynb
      
      c. Codon-level analysis - RiboSeq_CodonAnalysis.ipynb


