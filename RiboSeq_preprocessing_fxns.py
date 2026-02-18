import os
import re
import sys
import time
import json
import pysam
import subprocess
import multiprocessing
import matplotlib.pyplot as plt
import pandas as pd



##################### Call to paired .sh (bash) files for running appropriate command-line tools ########################

def runCutadapt(infile, cutfile, shortfile, longfile, cutlog, threePrimeAdapter, minLen, maxLen, acceptableError, numCores):
    '''
    Trims the provided 3' adapter from all reads (both full & partial occurences, plus all downstream nucleotides)
    Where resulting trimmed read is not in the specified length range, saved to appropriate long/short file
    
    Input:
        infile - full path for input data fastq or fastq.gz
        cutfile - full path for output (1) fastq.gz containing all passing reads
        shortfile - full path for output (2) fastq.gz containing all reads below minLen
        longfile - full path for output (3) fastq.gz containing all reads above maxLen
        cutlog - full path for output (4) txt file containing cutadapt report
        threePrimeAdapter - nucleotide sequence of 3' adapter to be trimmed
        minLen - minimum lenght of reads post adapter trimming (lower end of ribosome footprint size)
        maxLen - maximum lenght of reads post adapter trimming (upper end of ribosome footprint size)
        acceptableError - if N>1: total allowable errors in adapter matching; if 0≤N<1: acceptable ratio of # errors to adapter length 
        numCores - parallel processing capacity - must match scheduler
        
    Output:
        1. _cutadapt.fastq.gz - All reads with adapter trimmed and which pass the given length/quality filters. 
            Reads without adapters are kept until length processing.
        2. _cutshort.fastq.gz -all reads below the minLen after adapter trimming
        3. _cutlong.fastq.gz - all reads longer than maxLen after adapter trimming
        4. _cutlog.txt - cutadapt summary report
        
    cutadapt documentation @ https://cutadapt.readthedocs.io/en/stable/guide.html#basic-usage
    '''
    with open(cutlog, 'a') as log:
        log.write('\n ---cutadapt--- \n')
    subprocess.call(['bash', 'RiboSeq_runCutadapt.sh', infile, cutfile, shortfile, longfile, cutlog, threePrimeAdapter, str(minLen), str(maxLen), str(acceptableError), str(numCores)])


def runFastQC(adapterfile, outdir, infile):
    
    '''
    Input: 
        adapterfile: full path for text file containing relevant adapter names & sequences separated by tab
        outdir: folder to write FastQC results
        infile: full path for fastq or fastq.gz file to be analyzed
    Output: 
        html & .zip files containing FastQC analysis detailes
        html file can be opened in browser
    
    FastQC documentation @ https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/    
    '''
    
    subprocess.call(['bash', 'runFastQC.sh', adapterfile, outdir, infile])
    

def runBowtieBuild(reference_in, ebwt_base):
    '''
    build a Bowtie index from a set of DNA sequences
    Input: 
        reference_in - reference genome, in fasta format
        ebwt_base - name for output of the reference genome index
    Output: 
        outputs a set of 6 files with suffixes .1.ebwt, .2.ebwt, .3.ebwt, .4.ebwt, .rev.1.ebwt, and .rev.2.ebwt. 
        If the total length of all the input sequences is greater than about 4 billion, then the index files will end in ebwtl instead of ebwt 
        These files together constitute the index: they are all that is needed to align reads to that reference
    documentation: https://bowtie-bio.sourceforge.net/manual.shtml#the-bowtie-aligner
    '''
    subprocess.call(['bash', 'runBowtieBuild.sh', reference_in, ebwt_base])


def runBowtieAlign(fname, infile, filtering_outdir, genomeAlign_outdir, filtering_ebwt, genome_ebwt, bowtieLog, acceptableError, numCores):
    
    '''      
    updated 20240312 to run both SAM and BWT alignments
    bowtie takes an index and a set of reads as input and outputs a list of alignments 
    works best when aligning short reads to large genomes 
    
    Input:
        filtering/genome_ebwt: path to reference genome index or to index for tRNA/rRNA/etc filtering & removal
        infile: path to fastq file being aligned
        filtering/genomeAlign_outdir: path to alignment directories, separating filtering and final alignment results
        bowtieLog - name.extension for adding into log (not path)
    Output: 
        Aligns infile to ref filtering_ebwt for removal, then aligns filtered output to full genome_ebwt in SAM format
        Reads will also be listed in fastq files (_unalign.fastq, _align.fastq, _max.fastq), depending on alignment status with given conditions
        Note: these will be converted to fastq.gz
        logfile lists time and statistics data
    
    bowtie documentation: https://bowtie-bio.sourceforge.net/manual.shtml#the-bowtie-aligner      
    '''
    
    with open(bowtieLog, 'a') as log:
        log.write('\n ---BowtieAlign--- \n')
    subprocess.call(['bash', 'RiboSeq_runBowtieAlign.sh', fname, infile, filtering_outdir, genomeAlign_outdir, filtering_ebwt, genome_ebwt, bowtieLog, str(acceptableError), str(numCores)])




########### Read SAM files for WIG conversion & alignment QC #################

def SamToWig(fname, infile, outdir, gff_ID, minLen, maxLen):
    lenindex = range(minLen, maxLen+1)

    forward = {}
    reverse = {}
    fiveprime_fwd = {}
    fiveprime_rev = {}
    threeprime_fwd = {}
    threeprime_rev = {}
    
    samfile = pysam.AlignmentFile(infile, 'rb')
    for alignment in samfile.fetch():
        '''
        qName = alignment.query_name
        isforward = alignment.is_forward
        qSequence =alignment.query_sequence
        leftpos = alignment.reference_start #returns base 0 positioning
        NumMismatch = alignment.get_tag("XA")
        Mismatch = alignment.get_tag("MD")
        if x.get_tag("XM") !=0:
            print(x.get_tag("XA"), x.get_tag("MD")) 
        '''
        if alignment.is_mapped:
            
            seq = alignment.query_sequence #reverse complement if rev strand
            leftmost = alignment.reference_start + 1 #want 1-based offset into the forward strand (pysam.reference_start gives base 0)
            length = alignment.query_length
            rightmost = alignment.reference_end
            strand = alignment.is_forward #True if fwd strand, False if rev strand

            #list mismatch positions where MD=20T0G0 becomes ['20', '0', '0'], MD=0T0G20 becomes ['0', '0', '20']
            #GGCCAAGCAGCGTTGCCTTGNA
            #because MD lists position from last mutation, easier to count dist from start of read for 1st mismatch & dist from end of read for 2nd
            MD = alignment.get_tag('MD')
            mismatch_posList = re.split(r'\D+',MD)
            mismatchCount = len(mismatch_posList) - 1

            if mismatchCount ==1:
                mut_base = int(mismatch_posList[0])
            
                if mut_base == 0:
                    length -= 1
                    leftmost += 1   # move the 5' end of read (not yet accting for strand) 1 nt into fwd strand
                
                elif mut_base == length - 1:
                    length -= 1
                    rightmost -= 1
            
            elif mismatchCount == 2:
                mut_base_1 = int(mismatch_posList[0]) 
                mut_base_2 =  -(int(mismatch_posList[2])+1)
                
                if mut_base_1 == 0:
                    length -=1
                    leftmost += 1
    
                    if (mut_base_2 + length) == 0:
                        length -= 1
                        leftmost += 1
    
                if mut_base_2 == -1:
                    length -= 1
                    rightmost -= 1
    
                    if mut_base_1 == length-1:
                        length -= 1
                        rightmost -= 1

            
            if length in lenindex:
            #if (minLen<= length <= maxLen):
                if strand == True:
                    for x in range(leftmost, rightmost +1):
                        if x in forward:
                            forward[x] += 1.0/length
                        else:
                            forward[x] = 1.0/length
                            
                    if leftmost in fiveprime_fwd:
                        fiveprime_fwd[leftmost] += 1.0
                    else:
                        fiveprime_fwd[leftmost] = 1.0

                    if rightmost in threeprime_fwd:
                        threeprime_fwd[rightmost] += 1.0
                    else:
                        threeprime_fwd[rightmost] = 1.0

                elif strand == False:
                    for x in range(leftmost, rightmost + 1):
                        if x in reverse:
                            reverse[x] += 1.0 / length
                        else:
                            reverse[x] = 1.0 / length
                    if rightmost in fiveprime_rev:
                        fiveprime_rev[rightmost] += 1.0
                    else:
                        fiveprime_rev[rightmost] = 1.0
                    if leftmost in threeprime_rev:
                        threeprime_rev[leftmost] += 1.0
                    else:
                        threeprime_rev[leftmost] = 1.0

    samfile.close()
    
    #output 6 wig files
    writeoutwig(forward, outdir + fname + '_avg' + '_fwd_fromSam.wig', gff_ID)
    writeoutwig(reverse, outdir + fname + '_avg' + '_rev_fromSam.wig', gff_ID)

    writeoutwig(fiveprime_fwd, outdir + fname + '_5pr' + '_fwd_fromSam.wig', gff_ID)
    writeoutwig(fiveprime_rev, outdir + fname + '_5pr' + '_rev_fromSam.wig', gff_ID)

    writeoutwig(threeprime_fwd, outdir + fname + '_3pr' + '_fwd_fromSam.wig', gff_ID)
    writeoutwig(threeprime_rev, outdir + fname + '_3pr' + '_rev_fromSam.wig', gff_ID)


def readInfo_sam(fname, infile, json_dir, logfile, minLen, maxLen): 
    ### Set initial variables ###    
    
    print(fname, 'readInfo_sam')
    lenindex = range((minLen), (maxLen+1))

    lengthCounts = {length : 0 for length in lenindex}
    
    totalReads = 0

    '''
    ex: G_data = { length: [# reads of length with G in pos 0, # reads of length with G pos 1, ... up to maxLen],
                            # reads=0 for position > length
                    3: [0, 0, 0, ...],
                    4: [0, 0, 0, 0, ...],
                    ... }
                    
    --> each length has a list with positional composition of reads, starting at the 3' end
    '''
    
    G_data = {length : [0]*(maxLen+1) for length in lenindex}
    C_data = {length : [0]*(maxLen+1) for length in lenindex}
    A_data = {length : [0]*(maxLen+1) for length in lenindex}
    U_data = {length : [0]*(maxLen+1) for length in lenindex} 
    
    samfile = pysam.AlignmentFile(infile, 'rb')
    for alignment in samfile.fetch():
        '''
        qName = alignment.query_name
        isforward = alignment.is_forward
        qSequence =alignment.query_sequence
        leftpos = alignment.reference_start #returns base 0 positioning
        NumMismatch = alignment.get_tag("XA")
        Mismatch = alignment.get_tag("MD")
        if x.get_tag("XM") !=0:
            print(x.get_tag("XA"), x.get_tag("MD")) 
        '''
        if alignment.is_mapped:
            
            seq = alignment.query_sequence #reverse complement if rev strand
            leftmost = alignment.reference_start + 1 #want 1-based offset into the forward strand (pysam.reference_start gives base 0)
            length = alignment.query_length
            rightmost = alignment.reference_end
            strand = alignment.is_forward #True if fwd strand, False if rev strand

            read =alignment.get_forward_sequence()

            #list mismatch positions where MD=20T0G0 becomes ['20', '0', '0'], MD=0T0G20 becomes ['0', '0', '20']
            #GGCCAAGCAGCGTTGCCTTGNA
            #because MD lists position from last mutation, easier to count dist from start of read for 1st mismatch & dist from end of read for 2nd
            MD = alignment.get_tag('MD')
            mismatch_posList = re.split(r'\D+',MD)
            mismatchCount = len(mismatch_posList) - 1
  
            if mismatchCount ==1:
                mut_base = int(mismatch_posList[0])
            
                if mut_base == 0:
                    length -= 1
                    leftmost += 1   # move the 5' end of read (not yet accting for strand) 1 nt into fwd strand
                    
                    if strand ==True:
                        read = read[1:]
                    else:
                        read = read[:-1]
                
                elif mut_base == length - 1:
                    length -= 1
                    rightmost -= 1

                    if strand ==True:
                        read = read[:-1]
                    else: 
                        read = read[1:]

              
            
            elif mismatchCount == 2:
                mut_base_1 = int(mismatch_posList[0]) 
                mut_base_2 =  -(int(mismatch_posList[2])+1)
                
                
                if mut_base_1 == 0:
                    length -=1
                    leftmost += 1

                    if strand ==True:
                        read = read[1:]
                    else:
                        read = read[:-1]
                    
                    
    
                    if (mut_base_2 + length) == 0:
                        length -= 1
                        leftmost += 1

                        if strand ==True:
                            read = read[1:]
                        else:
                            read = read[:-1]

                    


        
                if mut_base_2 == -1:
                    length -= 1
                    rightmost -= 1

                    if strand ==True:
                        read = read[:-1]
                    else: 
                        read = read[1:]
    
                    if mut_base_1 == length-1:
                        length -= 1
                        rightmost -= 1

                        if strand ==True:
                            read = read[:-1]
                        else: 
                            read = read[1:]

        
            if ((minLen)<= length <= maxLen):

                totalReads += 1
                lengthCounts[length] += 1

            
                for position in range(0, length):

                    if read[position] == 'G':
                        G_data[length][position] += 1
                        
                    elif read[position] == 'A':
                        A_data[length][position] += 1
                        
                    elif read[position] == 'C':
                        C_data[length][position] += 1
                        
                    elif read[position] == 'T':
                        U_data[length][position] += 1    
    
    samfile.close()
    
    
    for length in lenindex:
        if lengthCounts[length] >= 1:
            tot = lengthCounts[length]   # tot= total number of nucleotides at any position for a given length
        else:
            tot = 1

        # transform absolute count to fraction of total in a given length of reads
        for position in range(0, maxLen):
            g = float(G_data[length][position]) / float(tot) *100
            c = float(C_data[length][position]) / float(tot) *100
            u = float(U_data[length][position]) / float(tot) *100
            a = float(A_data[length][position]) / float(tot) *100
    
             # replace absolute counts with fractions calculated above:
        
            G_data[length][position] = g
            C_data[length][position] = c
            U_data[length][position] = u
            A_data[length][position] = a


    nt_dict = {'G': G_data, 'C': C_data, 'U': U_data, 'A': A_data}

     # calculate total reads in library - sanity check for total reads
    
    sum_reads = sum(lengthCounts.values())
    sum_reads = float(sum_reads)

    #print(sum_reads)
    #print(totalReads, sum_reads == float(totalReads))
    
     

    lengthPercent = {length : 0 for length in lenindex}
    for length in lengthCounts.keys():
            lengthPercent[length] = lengthCounts[length]/totalReads *100
    
    ### Output ###
        
    json_dump(lengthCounts, json_dir + fname + '_readCounts_sam.json')
    json_dump(lengthPercent, json_dir + fname + '_readPercent_sam.json')
    json_dump(G_data, json_dir + fname + '_comp_G')
    json_dump(C_data, json_dir + fname + '_comp_C')
    json_dump(A_data, json_dir + fname + '_comp_A')
    json_dump(U_data, json_dir + fname + '_comp_U')

    #print(sorted(lengthCounts.items()))
    #print(sorted(lengthPercent.items()))

    with open(logfile,'a+') as readcomp_log:
        readcomp_log.write('\n'+ fname + ' read composition\n')
        readcomp_log.write('total reads: %s\n' % totalReads)
        readcomp_log.write('read distribution (absolute counts): %s\n' % sorted(lengthCounts.items()))
        readcomp_log.write('read distribution (%% of total reads): %s\n' % sorted(lengthPercent.items()))
        readcomp_log.write('G_comp: \n' + str(G_data))
        readcomp_log.write('\nC_comp: \n' + str(C_data))
        readcomp_log.write('\nA_comp: \n' + str(A_data))
        readcomp_log.write('\nU_comp: \n' + str(U_data))
                           

    
    
    return lengthCounts, lengthPercent, nt_dict





##################### Plot read length & nt composition ##############################

def plot_ReadLength(readQC_dict, outdir):
    
    for fname in readQC_dict.keys(): 

        readCount_dict = {int(k):int(v) for k,v in readQC_dict[fname]['lenCounts'].items()}
        readPercent_dict = {int(k):float(v) for k,v in readQC_dict[fname]['lenPercent'].items()} 
        
        readCounts = dict_to_df(readCount_dict, 'Read Length', 'Read Counts')
        readPercents = dict_to_df(readPercent_dict, 'Read Length', '% Total Reads')
        
        
        ax = readCounts.plot.bar()
        max_index, max_val = max(zip(readCounts.index, list(readCounts['Read Counts'])), key=lambda x: x[1])
        ax.bar(x=(max_index-min(readCount_dict.keys())), height=max_val, color='r', width=[p.get_width() for p in ax.patches][0])
        ax.set_title('{} Size Distribution, mode at {}'.format(fname, max_index))
        ax.set_xlabel('Read Length')
        ax.set_ylabel('Read Counts')
        ax.get_yaxis().get_major_formatter().set_scientific(False)
        plt.savefig( outdir + fname + '_readCounts.pdf', bbox_inches='tight')
        plt.close()
        
        
        ax = readPercents.plot.bar()
        max_index, max_val = max(zip(readPercents.index, list(readPercents['% Total Reads'])), key=lambda x: x[1])
        ax.bar(x=(max_index-min(readPercent_dict.keys())), height=max_val, color='r', width=[p.get_width() for p in ax.patches][0])
        ax.set_title('{} Size Distribution, mode at {}'.format(fname, max_index))
        ax.set_xlabel('Read Length')
        ax.set_ylabel('% of Total Reads')
        plt.savefig( outdir + fname + '_readPercents.pdf', bbox_inches='tight')
        plt.close()

    for fname in readQC_dict.keys(): 

        readCount_dict = {int(k):int(v) for k,v in readQC_dict[fname]['lenCounts'].items()}
        readPercent_dict = {int(k):float(v) for k,v in readQC_dict[fname]['lenPercent'].items()}
        
        readCounts = dict_to_df(readCount_dict, 'Read Length', 'Read Counts')
        readPercents = dict_to_df(readPercent_dict, 'Read Length', '% Total Reads')
        
        plt.plot(readPercents, label = fname)
        plt.title('Size Distribution')
        plt.xlabel("Read Length")
        plt.ylabel("Percent of Reads")
        lgd = ax.legend(fontsize=16, loc='center left', bbox_to_anchor=(1, 0.5))
        
    plt.savefig( outdir + 'compare_sizeDist_' + '.pdf', bbox_inches="tight")
    plt.close()

    

def plot_nucleotideComposition(readQC_dict, outdir, minLen, maxLen):
    #plot nucleotide composition of reads in dataset
    # nt_dict = {'G': G_data, 'C': C_data, ...} where G_data ={readlength: [G counts by position]}
    for fname in readQC_dict:
        
        nt_dict = readQC_dict[fname]['nt_dict']
        
        plot_num = 0
        plt.figure(figsize=(20,2.5))
        for nucleotide in ['G', 'C', 'A', 'U']:
            plot_num += 1
            print(fname + '_readcomp-' + nucleotide)
            
            data = heatmapdict_to_df(nt_dict[nucleotide], 'ReadLength', 'ReadPosition',('%s_Counts' % nucleotide))
            
            plt.subplot(1,4,plot_num)
            plt.style.context('white')
            
            plot = plt.imshow(data, cmap = "RdBu_r", vmin = 0, vmax=50)
            
            plt.colorbar()
            l, r = plt.xlim(); t, b = plt.ylim()
            plt.xticks(ticks = range(int(l), int(r)+1, 15), labels = range(0, maxLen+1, 15))
            print((int(t)-int(b)), (maxLen-minLen))
            plt.yticks(range(int(b), int(t)+1, 15) , range(minLen, maxLen+1, 15) )
           
            
            plt.xlabel('ReadPosition')
            plt.ylabel('ReadLength')
              
            plt.title(fname + ' - ' + nucleotide)
            
            
        plt.savefig(outdir + fname + '_read_composition.pdf', dpi=400, bbox_inches="tight") 
        plt.close()


####################### utility functions ########################

def writeoutwig(data, outFile, gff_ID):
    data_sort = {k: data[k] for k in sorted(data.keys())}
    outFile = open(outFile, 'w')
    outFile.write('track type=wiggle_0' + '\n')
    outFile.write('variableStep chrom='+ gff_ID + '\n')
    for y in list(data_sort.items()):
        base = str(y[0])
        read = str(y[1])
        outFile.write(base + '\t' + read + '\n')
    outFile.close()

def json_dump(data, resultsdir):
    if resultsdir[-5:] == '.json':
        with open(resultsdir , 'w') as f:
            json.dump(data, f)
    else:  
        with open(resultsdir + '.json', 'w') as f:
            json.dump(data, f)

def json_load(file):
    if file[-5:] == '.json':
        with open(file , 'r') as j:
            return json.load(j)
    else:
        with open(file + '.json', 'r') as j:
            return json.load(j)

def dict_to_df(dictionary, col_name, val_name):
    '''Make dictionary = {index: list_index(values)
    into pandas dataframe'''

    df = pd.DataFrame(columns = (col_name, val_name))
    col_index = sorted(dictionary.keys())
    collist   = []
    vallist   = []
    
    for column in col_index:
        collist.append(column)
        vallist.append(dictionary[column])      
        
    data = {val_name : vallist}
    df   = pd.DataFrame(data, index = collist)
   
    return df

def heatmapdict_to_df(dict, col_name, row_name, val_name):
    df = pd.DataFrame(columns = (col_name, row_name, val_name))
    col_index = sorted(dict.keys())
    row_index = [x for x in range(0, len(dict[col_index[1]]))]
    
    collist   = []
    rowlist   = []
    vallist   = []

    for column in col_index:
        for row in row_index:
            collist.append(column)
            rowlist.append(row)
            if dict[column][row] > 0:
                vallist.append(dict[column][row])
            else:
                vallist.append(None)
            

    data = {col_name : collist, 
            row_name : rowlist,
            val_name : vallist
           }
    df = pd.DataFrame(data)
    df = df.pivot(index=col_name, columns = row_name, values = val_name)
    return df
        