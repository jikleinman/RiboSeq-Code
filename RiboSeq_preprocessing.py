#!/usr/bin/env python
# coding: utf-8

# # RiboSeq Preprocessing
# Performs trimming and alignment functions on both/either RiboSeq & RNAseq reads
# 
#     - input: demultiplexed sequencing fastq.gz files, stored in a subdirectory seq_gz
#     - output: QC, trimmed files, alignment files, WIG files (used for data analysis steps)
# 
# convert to .py before running from wynton scheduler!
# 
# 8 samples takes ~1.5-2.5h to run on wynton server
# 
#     from terminal: 
#     
#     # create python script
#     > jupyter nbconvert --to python RiboSeq_preprocessing.ipynb
#     
#     # run this script from commandline directly or submit to remote server
#     # Not recommended to run directly in jupyter notebook
#     
#     # fill out RiboSeq_preprocessing_scheduler.sh and submit by:
#     > qsub RiboSeq_preprocessing_scheduler.sh
# 
#     # or run directly in terminal (conda environment with dependencies must be active!!)
#     > python RiboSeq_preprocessing.py filename1 filename2...
#     

# #### version info

# version 5 = 
# 
# 2024/09/18
# 
#     added QC with WIG conversion trimming - nucleotide composition plots & size distribution/size distribution comparison plots
# 
# 2024/09/11 JIK 
#     
#     added QC step after alignment,  added specific rRNA/tRNA/size control filtering step during bowtie alignment, removed bowtie BWT outputs to save processing time/power, 
# 
#     
# version 4 = 2024/07/05 JIK

# ## 1. Importing useful modules & defining necessary functions
# 
# Scroll down to the box beginning with "if __name__ == '__main__':" for user inputs
# 

# In[ ]:


# RiboSeq Preprocessing 
import os
import re
import sys
import time
import pysam
import subprocess
import multiprocessing

import RiboSeq_preprocessing_fxns as preprocessing

from importlib import reload
#reload(preprocessing)


# ## 2. User-defined inputs
# 
# Check these variables closely!! 
# * Adapters, dapterfile, reference genome, and reference filter sequences need to match your experimental design
#   
# * Footprint lengths can match experimental design or be more stringent
# 
# * dir variables tell the code where to find input info and where to spit out the output
# 
#    * code will throw an error if these don't match your directory architecture! <p>
# 
# * acceptable error represents quality filtering
# 
#    * acceptableError = 2 means that a maximum of 2 mismatches are allowed in the adapter region for trimming, or else the read is excluded AND that reads with >2 mismatches during alignment to genome will be excluded <p>
# 
# * numCores is used for parallel processing - allowing you to analyze multiple data files at once. This variable MUST match the one in RiboSeq_preprocessing_scheduler.sh
#      

# In[ ]:


if __name__ == '__main__':
    
    #file_names = sys.argv[1:]
    file_names = ['JKDF06_Rn1_10xDZD_S71_L003_R1_001']
    
    ###################################################################################################################################
    ### check that all information here matches experimental details, folder layout, & scheduler ######################################
    ###################################################################################################################################
    
    # overarching directory - code will output into a subfolder of this location
    maindir = '/wynton/home/fujimori/jkleinman/JK14075/' #'/wynton/home/fujimori/jkleinman/JKDF05/'
    
    # directory where raw data is saved - can be subdirectory of maindir or elsewhere
    raw_data_dir = maindir + 'seq_gz/'
    
    # directory where scripts and adapter files can be found - run code from within this directory
    # code is generally stored in a folder separate from data, but can be either way
    codedir = '/wynton/home/fujimori/jkleinman/RiboSeq_code/'
    
    #text file containing relevant adapters used in the riboseq protocol (name + adapter, tab separated)
    #used as a confirmation of cutadapt success in QC
    fastqc_adapterfile = codedir + 'RiboSeq_adapters.txt'
    
    # directory for bowtie alignment indexes (if this folder doesn't already exist, code will build it)
    BowtieIndex_dir = codedir + 'indexes/'
    
    # directory for output of data processing & analysis - code builds this folder for you, dated by when the code is run
    preprocessing_outdir = maindir + time.strftime("%Y%m%d") + '_RiboSeq-preprocessing_%s/'
    
    # 3' adapter removed during cutadapt
    # NEB universal miRNA cloning primer = 'CTGTAGGCACCATCAAT'
    threePrimeAdapter = 'CTGTAGGCACCATCAAT'
    
    # minimum length of ribosome footprint
    minLen = 20
    
    # maximum length of ribosome footprint
    maxLen = 100
    
    # allowable number of errors per read in adaptor matching during cutadapt & bowtie alignment
    acceptableError = 2
    
    # name and file location of reference genome for alignment (make sure you have at least the .fasta file in here, but you'll also want the .gff later)
    # Best practice is to rename your reference genome files & place them in a genomefiles subfolder with the same name before running the code
    # I prefer the format [genome].[versionNumber] for consistency 
    reference_genome_ID = 'CP009273.1'
    # this should be the folder where you've saved the fasta & gff files for your genome
    # eg. folder RS_code/genomefiles/CP009723.1/ (which is what the below line says) must contain file CP009723.1.fasta
    reference_genome_fasta = codedir + 'genomefiles/CP009273.1/CP009273.1.fasta'

    #  before running this code, create a fake genome file containing  all of the rRNA & tRNA sequences in your genome. extra linker DNA not necessary.
    # if you have other sequences you wish to filter, simply add them into the same fasta file before running code
    # This will be used as a filter during alignment; give the ID a descriptive name since this will show on the log as what is being filtered out.
    # see 'CP009273.1_tRNArRNA' files for reference
    reference_filtering_ID = 'CP009273.1_tRNArRNAtmRNAladderillumina'
    reference_filtering_fasta = codedir + '/genome_files/CP009273.1/CP009273.1_tRNArRNAtmRNAladderillumina.fa'

    # number of cores to be used for parallel processing. Must match value in RiboSeq_preprocessing_scheduler.sh
    numCores = 40

    '''
    # Current quality analysis is based on # errors in matching & alignment, may want to add a PHRED filter back in at some point
    # Avg PHRED quality needed for read to pass preprocessing
    avgQuality = 30
    '''

    ###################################################################################################################################
    ###################################################################################################################################
    
    print('input: ', file_names)
    


# ## 3. Preprocessing
# 
# This shouldn't generally require any user input, but this is all the actual data manipulation to remove adapters & align to genome so it's worth understanding what's going on

# In[ ]:


# timing
    starttime = time.time()


# In[ ]:


# create outdir for total preprocessing; all other output folders will be nested in this directory
    i=1
    if not os.path.exists(preprocessing_outdir[:-4] + '/'):
        preprocessing_outdir = preprocessing_outdir[:-4] + '/'
    else:
        while os.path.exists(preprocessing_outdir % str(i)):
            i+=1
        preprocessing_outdir = preprocessing_outdir % str(i)
    
    print(preprocessing_outdir)
    os.mkdir(preprocessing_outdir)


# In[ ]:


# check input file type (fastq vs. fastq.gz)
# create a dictionary of all raw files where {fname: full path of file}
    
    rawfiles = {}

    for fname in file_names:

        infile = raw_data_dir+fname+'.fastq'
        
        if os.path.exists(infile):
            print(os.path.split(infile)[1])
            rawfiles[fname] = infile
        elif os.path.exists(infile + '.gz'):
            infile = infile+'.gz'
            print(os.path.split(infile)[1])
            rawfiles[fname] = infile
        else:
            print(infile +' path does not exist')


    for fname in rawfiles:
        print(fname, rawfiles[fname])


# In[ ]:


# create a log to summarize output from each step for a given data file

    # check for & create log directory
    logdir = preprocessing_outdir + 'logs/'
    if not os.path.exists(logdir):
        os.mkdir(logdir)

    # create a dictionary of logfiles for downstream referencing where {fname: full path of logfile}
    logfiles ={}

    for fname in rawfiles:
        logfiles[fname] = logdir + fname + '_log.txt'
        with open(logfiles[fname], 'a') as log:
            log.write('Preprocessing log: %s \n' % fname)
        


# In[ ]:


# run fastQC on raw files - outputs an html file for visual analysis by the user

    # check for output directory
    rawQC_outdir = preprocessing_outdir + 'rawQC/'
    if not os.path.exists(rawQC_outdir):
        os.mkdir(rawQC_outdir)
        
    
    # prepare arguments for all input files to be passed to fastQC    
    t1 = time.time()
    print('\nfastqc - raw files')
    rawQC_pool = multiprocessing.Pool(processes=numCores)
    rawQC_args = []
    
    for fname in rawfiles:
        
        infile = rawfiles[fname]
        
        rawQC_args.append((fastqc_adapterfile, rawQC_outdir, infile))
    
    # perform fastQC (runFastQC, calls runFastQC.sh)
    rawQC_results = rawQC_pool.starmap(preprocessing.runFastQC, rawQC_args)
    rawQC_pool.close()
    rawQC_pool.join()
    
    t2 = time.time()
    print('\nraw fastQC completed in', time.strftime("%H:%M:%S", time.gmtime(t2-t1)), 'hh:mm:ss\n')  
     


# In[ ]:


# run cutadapt on raw files (trim 3' adapters)

    # check for output directory
    cutadapt_outdir = preprocessing_outdir + 'cutadapt/'
    if not os.path.exists(cutadapt_outdir):
        os.mkdir(cutadapt_outdir)
    
    # prepare a dictionary {fname: full path of file} for cutadapt output fastq.gz containing all reads which passed QC
    cutfiles = {}
    
    
    # prepare arguments for all input files to be passed to cutadapt
    t1 = time.time()
    print('cutadapt')
    cut_pool = multiprocessing.Pool(processes=numCores)
    cut_args = []
    for fname in rawfiles:
        
        infile = rawfiles[fname]
        
        cutfile = cutadapt_outdir +fname + '_cutadapt.fastq.gz'
        shortfile = cutadapt_outdir +fname +'_cutshort.fastq.gz'
        longfile = cutadapt_outdir +fname +'_cutlong.fastq.gz'
        log = logfiles[fname]
        
        cutfiles[fname]=cutfile
        
        cut_args.append((infile, cutfile, shortfile, longfile, log, threePrimeAdapter, minLen, maxLen, acceptableError, numCores))
        
    # perform cutadapt (RiboSeq_runCutadapt, calls RiboSeq_runCutadapt.sh)
    cut_results = cut_pool.starmap(preprocessing.runCutadapt, cut_args)
    cut_pool.close()
    cut_pool.join()
    
    
    t2 = time.time()
    print('total cutadapt completed in ', time.strftime("%H:%M:%S", time.gmtime(t2-t1)), 'hh:mm:ss\n')
    


# In[ ]:


# run fastQC on adapter-trimmed files
    
    # check for output directory
    cutQC_outdir = preprocessing_outdir + 'cutQC/'
    if not os.path.exists(cutQC_outdir):
        os.mkdir(cutQC_outdir)
        
    # prepare arguments for all input files to be passed to fastQC    
    t1 = time.time()
    print('fastqc - cut files')
    cutQC_pool = multiprocessing.Pool(processes=numCores)
    cutQC_args = []
    
    for fname in cutfiles:
        
        infile = cutfiles[fname]
        
        cutQC_args.append((fastqc_adapterfile, cutQC_outdir, infile))
    
    # perform fastQC (runFastQC, calls runFastQC.sh)
    cutQC_results = cutQC_pool.starmap(preprocessing.runFastQC, cutQC_args)
    cutQC_pool.close()
    cutQC_pool.join()

    t2 = time.time()
    print('cutadapt fastQC completed in ', time.strftime("%H:%M:%S", time.gmtime(t2-t1)), 'hh:mm:ss\n')
     


# In[ ]:


# Check for & create bowtie indices

# if .ebwt index files already exist, those will be used instead of rebuilding the index
# to rebuild .ebwt files, run the commented out part of this block

    # check for BowtieIndex directory
    if not os.path.exists(BowtieIndex_dir):
        os.mkdir(BowtieIndex_dir)


    # build a new bowtie index - 1st time only 
    # checks whether .ebwt files already exist (no user input needed unless you want to rebuild them for some reason)

    if not os.path.exists(BowtieIndex_dir + reference_genome_ID + '.1.ebwt'):
        
        # run bowtieBuild
        preprocessing.runBowtieBuild(reference_genome_fasta, reference_genome_ID)
       
    # repeat the process to get indexes for the tRNA/rRNA filtering fasta
    if not os.path.exists(BowtieIndex_dir + reference_filtering_ID + '.1.ebwt'):
    
        # run bowtieBuild
        preprocessing.runBowtieBuild(reference_filtering_fasta, reference_filtering_ID)
    
    # move the newly created indices to the proper location
    for fname in os.listdir(codedir):
        if '.ebwt' in fname:
            os.rename(codedir+fname, BowtieIndex_dir+fname)

    
    # run this chunk to rebuild (rewrite!) an existing genome index
    '''
    preprocessing.runBowtieBuild(reference_genome_fasta, reference_genome_ID)
    for fname in os.listdir(codedir):
        if '.ebwt' in fname:
            os.rename(codedir+fname, BowtieIndex_dir+fname)
    '''
    # run this chunk to rebuild (rewrite!) an existing filtering index
    '''
    preprocessing.runBowtieBuild(reference_filtering_fasta, reference_filtering_ID)
    for fname in os.listdir(codedir):
        if '.ebwt' in fname:
            os.rename(codedir+fname, BowtieIndex_dir+fname)    
    '''


# In[ ]:


# run Bowtie Alignment on adapter-trimmed files (20240911 sam only since bwt and sam have been consistent)

# note: SAM file may be longer than BWT due to issue where some unaligned reads are placed in SAM file
# these reads have FLAG = 4 which indicates that they are not mapped to genome
# & will thus be filtered out by pysam.is_mapped during conversion to WIG

   
    # check for output directory
    bowtieAlign_outdir = preprocessing_outdir + 'bowtieAlign/'
    filtering_outdir = bowtieAlign_outdir + 'filter/'
    genomeAlign_outdir = bowtieAlign_outdir + 'genome_align/'
    
    
    if not os.path.exists(genomeAlign_outdir):
        os.mkdir(bowtieAlign_outdir)
        os.mkdir(filtering_outdir)
        os.mkdir(genomeAlign_outdir)
        
    # bowtie alignment

    genomeSamFiles = {}
    
    
    # prepare arguments for all input files to be passed to bowtie
    t1 = time.time()
    print('bowtie alignment')
    bowtie_pool = multiprocessing.Pool(processes=numCores)
    bowtie_args = []
    
    genome_ebwt = BowtieIndex_dir + reference_genome_ID
    filtering_ebwt = BowtieIndex_dir + reference_filtering_ID
    
    for fname in cutfiles:
        
        infile = cutfiles[fname]
        
        # will output 5 files: fname_align.sam, fname_unalign.fastq.gz, fname_max.fastq.gz, fname_align.fastq.gz, fname_bowtieLog.txt
        # .sam are the alignment file extensions to be used in downstream steps
        
        bowtieLog = logfiles[fname]
        alignSam = fname + '_filtered_genome_match.sam'
        genomeSamFiles[fname] = genomeAlign_outdir + alignSam
        
        bowtie_args.append((fname, infile, filtering_outdir, genomeAlign_outdir, filtering_ebwt, genome_ebwt, bowtieLog, acceptableError, numCores))
    
    # perform bowtie alignment (RiboSeq_runBowtieAlign, calls RiboSeq_runBowtieAlign.sh)
    bowtie_results = bowtie_pool.starmap(preprocessing.runBowtieAlign, bowtie_args)
    bowtie_pool.close()
    bowtie_pool.join()

    t2 = time.time()
    print('bowtie alignment completed in ', time.strftime("%H:%M:%S", time.gmtime(t2-t1)), 'hh:mm:ss\n') 

    '''for fname in cutfiles:
        infile = cutfiles[fname]
        alignSam = bowtieAlign_outdir + fname + '_align.sam'
        alignedfiles[fname] = alignSam'''
        


# In[ ]:


# run fastQC on bowtie aligned files
    
    # check for output directory
    bowtieQC_outdir = preprocessing_outdir + 'bowtieQC/'
    if not os.path.exists(bowtieQC_outdir):
        os.mkdir(bowtieQC_outdir)
        os.mkdir(bowtieQC_outdir + 'filtered/')
        os.mkdir(bowtieQC_outdir + 'genome_aligned/')
        
        
    # prepare arguments for all input files to be passed to fastQC    
    t1 = time.time()
    print('fastqc - alignment fastq files')
    bowtieQC_pool = multiprocessing.Pool(processes=numCores)
    bowtieQC_args = []
    
    for fname in genomeSamFiles:

        filtFastq = filtering_outdir + fname + '_filtering_removed.fastq.gz'
        genomeFastq = genomeAlign_outdir + fname + '_filtered_genome_match.fastq.gz'
        
        bowtieQC_args.append((fastqc_adapterfile, bowtieQC_outdir + 'filtered/', filtFastq))
        bowtieQC_args.append((fastqc_adapterfile, bowtieQC_outdir + 'genome_aligned/', genomeFastq))
    
    # perform fastQC (runFastQC, calls runFastQC.sh)
    bowtieQC_results = bowtieQC_pool.starmap(preprocessing.runFastQC, bowtieQC_args)
    bowtieQC_pool.close()
    bowtieQC_pool.join()

    t2 = time.time()
    print('bowtie fastQC completed in ', time.strftime("%H:%M:%S", time.gmtime(t2-t1)), 'hh:mm:ss\n')
     


# In[ ]:


#convert SAM to WIG

    print('WIG conversion -- minLen =%s' %minLen, ' maxLen=%s' % maxLen, '\nUsing minLen as cutoff after end mismatch trimming\n')
    #check for output directory
   
    wig_outdir = preprocessing_outdir + 'Wig/'
    
    if not os.path.exists(wig_outdir):
        os.mkdir(wig_outdir)
        
        
    t1 = time.time()
    samWig_pool = multiprocessing.Pool(processes=numCores)
    samWig_args = []
    for fname in genomeSamFiles:
        infile = genomeSamFiles[fname]
        with open(logfiles[fname], 'a') as log:
            log.write('\n ---WIG conversion - minLen=%s, maxLen=%s --- \n' %(minLen,maxLen))
            log.write('Using minLen as cutoff after end mismatch trimming\n')
        
        samWig_args.append((fname, infile, wig_outdir, reference_genome_ID, minLen, maxLen))
        
    wig_results = samWig_pool.starmap(preprocessing.SamToWig, samWig_args)
    samWig_pool.close()
    samWig_pool.join()
    
    t2 = time.time()
    print('Sam to Wig completed in ', time.strftime("%H:%M:%S", time.gmtime(t2-t1)), 'hh:mm:ss') 



# In[ ]:


# run custom QC on SAM files using same trimming as for WIG conversion
# output read length info (counts & percents) & nt composition by length (heat map of % G/C/A/U at each position in a read, by read length)

    json_dir = bowtieQC_outdir + 'genome_aligned/QC_json/'
    if not os.path.exists(json_dir):
        os.mkdir(json_dir)

    samQCdict = {}

    t1 = time.time()
    samQC_pool = multiprocessing.Pool(processes=numCores)
    samQC_args = []
    
    for fname in genomeSamFiles:
        infile = genomeSamFiles[fname]
        logfile = logfiles[fname]
        
        samQC_args.append((fname, infile, json_dir, logfile, minLen, maxLen))

        
    QC_results = samQC_pool.starmap(preprocessing.readInfo_sam, samQC_args)
    samQC_pool.close()
    samQC_pool.join()

    for fname in genomeSamFiles:
        
        samQCdict[fname] = {'lenCounts': preprocessing.json_load(json_dir + fname + '_readCounts_sam.json'), 
                            'lenPercent': preprocessing.json_load(json_dir + fname + '_readPercent_sam.json'), 
                            'nt_dict': {'G': preprocessing.json_load(json_dir + fname + '_comp_G'), 
                                       'C': preprocessing.json_load(json_dir + fname + '_comp_C'), 
                                       'A': preprocessing.json_load(json_dir + fname + '_comp_A'),
                                       'U': preprocessing.json_load(json_dir + fname + '_comp_U')
                                      }}
        
    preprocessing.plot_ReadLength(samQCdict, bowtieQC_outdir + 'genome_aligned/')
    preprocessing.plot_nucleotideComposition(samQCdict, bowtieQC_outdir + 'genome_aligned/', minLen, maxLen)
    
    t2 = time.time()
    print('Sam QC completed in ', time.strftime("%H:%M:%S", time.gmtime(t2-t1)), 'hh:mm:ss') 


# In[ ]:


# timing 
    endtime = time.time()
    print('\nPreprocessing complete on %d samples in ' %(len(file_names)) + time.strftime("%H:%M:%S", time.gmtime(endtime-starttime)), 'hh:mm:ss\n\n')


# In[ ]:


# 3. Data analysis!
    print('\nPreprocessing is now complete')
    print('\n* Use the FastQC files to look at quality in the actual reads with vs without adapters')
    print('\n* Check genome aligned FastQC, size distribution, and nucleotide composition plots for each sample')
    print('\n\t* FastQC -- check for overrepresented sequences (may need to add to filter file)')
    print('\n\t* nucleotide composition plots & size distribution comparison -- determine if any read lengths need to be removed from analysis')
    
    print('\nFrom here, move to the data analysis docs: ###########')


# ## 4. Data analysis!
# 
# Preprocessing is now complete 
# * Use the FastQC files to look at quality in the actual reads with vs without adapters
# * Check genome aligned FastQC, size distribution, and nucleotide composition plots for each sample
#     * FastQC -- check for overrepresented sequences (may need to add to filter file)
#     * nucleotide composition plots & size distribution comparison -- determine if any read lengths need to be removed from analysis
# 
# From here, move to the data analysis docs: ###########

# In[ ]:





# ## 5. Misc useful info 

# In[ ]:


#structure of BWT alignment file
    """
    structure of BWT file = 8+ columns of information.
        0. The read name, which is completely unique.
        1. The barcode info
        2. The sequence strand
        3. The genome to which the alignment was made.
        4. 0-based offset into the forward reference strand where leftmost character of the alignment occurs. (The position 5' of the most 5' base of the forward strand read or the position 3' of the 3' end of the reverse strand read).
        5. For forward strand, the read.  For reverse strand, the reverse complement of the read.
        6. Read quality score.
        7. Number of alternative aligments.
        8. Mismatches.

    Unlike SAM files, BWT files do not have a header
    IGV uses base 1 ; bwt base 0 ; SAM base 1
    """


# In[ ]:


#structure of SAM alignment file
    '''
    structure of SAM file = 11+ columns of information
    https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
    https://en.wikipedia.org/wiki/SAM_(file_format)

    1 	QNAME 	String 	Query template NAME (read name)
    2 	FLAG 	Int 	bitwise FLAG (flags during alignment)
    3 	RNAME 	String 	References sequence NAME (reference genome for alignment)
    4 	POS 	Int 	1- based leftmost mapping POSition 
            (1-based offset into the forward reference strand where leftmost character of the alignment occurs
            (for rev strand this will refer to 3' end, for fwd strand refer to 5' end))
    5 	MAPQ 	Int 	MAPping Quality (−10 log10 Pr{mapping position is wrong}, rounded to the nearest integer. 
            A value 255 indicates that the mapping quality is not available.)
    6 	CIGAR 	String 	CIGAR string (CIGAR string representation of alignment, eg 26M = 26 matches/mismatches)
    7 	RNEXT 	String 	Ref. name of the mate/next read (--only relevant for PE)
    8 	PNEXT 	Int 	Position of the mate/next read (-- only relevant for PE)
    9 	TLEN 	Int 	observed Template LENgth (--only relevant for PE)
    10 	SEQ 	String 	segment SEQuence (reverse-complemented if aligned to the reverse strand)
    11 	QUAL 	String 	ASCII of Phred-scaled base QUALity+33 (encoded by sequencer)
    12+ optional fields
        XM:i:<N> The number of mismatches in the alignment. Only present if SAM record is for an aligned read
        XA:i:<N> Aligned read belongs to stratum <N>. 
            In the -v alignment mode, an alignment’s stratum is defined as the total number of mismatches in the entire alignment
        MD:Z:<S> A string representation of the mismatched reference bases in the alignment. See SAM Tags format specification for details. Only present if SAM record is for an aligned read.
        NM:i:<N> The edit distance; that is, the minimal number of one-nucleotide edits (substitutions, insertions and deletions) needed to transform the read string into the reference string. Only present if SAM record is for an aligned read.

    bowtie flags= 
    1:The read is one of a pair
    2:The alignment is one end of a proper paired-end alignment
    4:The read has no reported alignments
    8:The read is one of a pair and has no reported alignments
    16:The alignment is to the reverse reference strand
    32:The other mate in the paired-end alignment is aligned to the reverse reference strand
    64:The read is mate 1 in a pair
    128:The read is mate 2 in a pair
    an unpaired read that aligns to the reverse reference strand will have flag 16. 
    A paired-end read that aligns and is the first mate in the pair will have flag 83 (= 64 + 16 + 2 + 1).

    CIGAR representation https://en.wikipedia.org/wiki/Sequence_alignment#Representations
    '''


# In[ ]:




