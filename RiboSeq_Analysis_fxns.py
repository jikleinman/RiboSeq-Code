
from Bio import SeqIO, Seq
import numpy as np
from collections import Counter
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import copy
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
import statsmodels.stats
import csv


aaList = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q',
          'C', 'G', 'P', 'A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V', '*']

gencode = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
             }

gencode = {k: v for k, v in sorted(gencode.items(), key=lambda x: x[1])}

aaTable = {}
for codon, aa in gencode.items():
    if aa not in aaTable:
        aaTable[aa] = []
    aaTable[aa].append(codon)

aaTable = {k: v for k, v in sorted(aaTable.items(), key=lambda x: x[0])}

################## functions ####################

def parse_genome(user_inputs):
    '''
    input: file paths for genome .gb and .fasta where .fasta contains the full genome sequence & .gb contains all file annotations
    output: dict of 3 dictionaries containing sequence, strand, and location info for non-pseudogene CDS only
    sequence dict contains full sequence with UTRs:                     full_sequence_dict = { featureName : [...ATG...TAG...] } 
    strand dict:                                                        strand_dict        = { featureName : 1, featureName2 : -1, ... }
    location dict contains genome positions of gene beginning & end:    location_dict      = { featureName : (start, end) } (not including UTR)
    '''
    genome_genbank = user_inputs['genome_genbank']
    genome_fasta = user_inputs['genome_fasta'] 
    utr_length_to_include = user_inputs['utr_length_to_include']

    genome_dict = {}
    full_sequence_dict = {}
    strand_dict = {}
    location_dict = {}
    #utr_sequence_dict = {}
    
    genome_seq = list(SeqIO.parse(genome_fasta, 'fasta'))[0]
    genbank_seq = list(SeqIO.parse(genome_genbank, 'genbank'))[0]
    
    for feature in genbank_seq.features:
        
        ### exclude non-coding features and pseudogenes from analysis
        if feature.type != 'CDS':
            continue
        elif 'pseudo' in feature.qualifiers:
            continue
    
        #### CDS features will be named by locus_geneName
        # ex: for CP009273.1, hns == BW25113_1237_hns
        name = feature.qualifiers['locus_tag'][0] + '_' + feature.qualifiers['gene'][0]
        start = feature.location.start
        end = feature.location.end

        ### extract sequences from genome_seq, including UTRs 
        # if utr_length_to_include =50, then gene sequence begins 50nt before start codon & ends 50 nt after stop codon
        ### positive strand genes
        if feature.location.strand == 1:
            seq = str(genome_seq.seq[start-utr_length_to_include:end+utr_length_to_include])
        ### negative strand genes (need to be reverse complemented)
        elif feature.location.strand == -1:
            seq = str(genome_seq.seq[start-utr_length_to_include:end+utr_length_to_include].reverse_complement())
        
        ### Pair gene sequence with feature name in dict
        full_sequence_dict[name] = seq

        '''
        ### Create a separate dictionary for 5'UTR in case it may be helpful
        utr_seq = seq[:utr_length_to_include]
        utr_sequence_dict[name] = utr_seq
        '''
        ### Store strand and locations as their own dictionaries
        strand_dict[name] = feature.location.strand
        location_dict[name] = (start, end)

    genome_dict = {'full_sequence_dict': full_sequence_dict, 'strand_dict': strand_dict, 'location_dict': location_dict, 
                   'total_genome_len': len(str(genome_seq.seq))}
    return genome_dict 

def read_in_wig_addOffset(sample_files_dict, genome_dict, user_inputs):

    ''' 
    input: sample_files_dict = { sample_name : [ fwd file_path, rev file_path ], sample2: [fwd, rev], ... }
            genome_dict = {'full_sequence_dict': {gene_name: [...ATG...TAG...], gene2: [],...} , 
                           'strand_dict':        {gene_name: 1, gene2: -1, ...}
                           'location_dict':      {gene_name: (start, end), gene2: (start, end), ...} }
    
     output: feature_dict_meta = { sample_name: {gene1:[reads by position], gene2: [], ...}
     '''
    
    utr_length_to_include = user_inputs['utr_length_to_include']
    total_genome_len = genome_dict['total_genome_len']

    feature_dict_meta = {}

    full_sequence_dict = genome_dict['full_sequence_dict']
    location_dict = genome_dict['location_dict']
    strand_dict = genome_dict['strand_dict']
    
    for sample_name in sample_files_dict:
        if len(user_inputs['rnaSeq_pairs'].keys()) != 0 and sample_name in user_inputs['rnaSeq_pairs'].values():
            mapping_offset = 0
        else:
            mapping_offset = user_inputs['mapping_offset']

        fwd, rev = sample_files_dict[sample_name]
        print('##### {} - offset {}'.format(sample_name, mapping_offset))
        feature_dict_meta[sample_name] = {}
        
        ### Load all read info from both fwd & rev .wig files, shifted by mapping_offset value
        # converts wig files to 2 dictionaries where {position_in_genome(shifted by mapping_offset): counts} for each of fwd & reverse strands
        fwd_dict = {}
        rev_dict = {}

        ##### forward #####
        
        with open(fwd) as infile:
            for line in enumerate(infile):
                if line[0] > 1:     ###Ignore the first TWO LINES = HEADER of the wig file
                    split_line = line[1].split('\t')
                    # negative mapping offset --> adj_pos can be negative or >= 0 and <total_len_genome
                    # positive mapping offset --> adj_pos can be > total_genome_len or >= 0 and <=total_len_genome
                    adjusted_position = int(split_line[0]) + mapping_offset                                        # Note: mapping offset addition
                    readCount = float(split_line[1])
                    # wig is sorted by position, so read in low to high
                    # circular genome --> negative pos map to end of genome, pos > genome length map to start of genome
                    if adjusted_position < 0:                         # negative mapping offset only
                        fwd_dict[total_genome_len + adjusted_position] = readCount
                    elif (0 <= adjusted_position <= total_genome_len):
                        if adjusted_position in fwd_dict.keys():
                            fwd_dict[adjusted_position] += readCount
                        else:
                            fwd_dict[adjusted_position] = readCount
                    elif adjusted_position > total_genome_len:         # positive mapping offset only
                        if (adjusted_position - total_genome_len) in fwd_dict.keys():
                            fwd_dict[adjusted_position - total_genome_len] += readCount
                        else:
                            fwd_dict[adjusted_position - total_genome_len] = readCount
                            
        print('Done with fwd')

        ##### reverse #####

        with open(rev) as infile:
            for line in enumerate(infile):
                if line[0] > 1:     ###Ignore the first TWO LINES = HEADER  of the file
                    split_line = line[1].split('\t')
                    # positive mapping offset --> adj_pos can be negative or >= 0 and <total_len_genome
                    # negative mapping offset --> adj_pos can be > total_genome_len or >= 0 and <=total_len_genome
                    adjusted_position = int(split_line[0]) - mapping_offset                                        # Note: mapping offset subtraction
                    readCount = float(split_line[1])
                    # wig is sorted by position, so read in low to high
                    # circular genome --> negative pos map to end of genome, pos > genome length map to start of genome
                    if adjusted_position < 0:                         # positive mapping offset only
                        rev_dict[total_genome_len + adjusted_position] = readCount
                    elif (0 <= adjusted_position <= total_genome_len):
                        if adjusted_position in rev_dict.keys():
                            rev_dict[adjusted_position] += readCount
                        else:
                            rev_dict[adjusted_position] = readCount
                    elif adjusted_position > total_genome_len:         # negative mapping offset only
                        if (adjusted_position - total_genome_len) in rev_dict.keys():
                            rev_dict[adjusted_position - total_genome_len] += readCount
                        else:
                            rev_dict[adjusted_position - total_genome_len] = readCount
                    
        print('Done with rev')
    
        ### get read info across genes defined in full_sequence_dict (includes UTRs)
        for gene_name in full_sequence_dict.keys(): 
            
            #### Start with positive strand genes 
            if strand_dict[gene_name] == 1:
                
                ### Get positions for gene start/end WITH UTR (-utr#, +utr#)
                pos = (location_dict[gene_name][0]-utr_length_to_include, location_dict[gene_name][1]+utr_length_to_include)
                if pos[1] < pos[0]:
                    print('found a bug: {} -- start-utr: {} end+utr: {} strand: {}'.format(gene_name, pos[0], pos[1], strand_dict[gene_name]))
                    continue
                    
                ### Create a list of read values across the length of gene with UTRs; 0 if there are no reads at that position
                sequencing = []
                for i in range(pos[0], pos[1]):
                    try:
                        sequencing.append(fwd_dict[i])
                    except KeyError:
                        sequencing.append(0)
                        
                feature_dict_meta[sample_name][gene_name] = sequencing
                
            #### repeat for negative strand genes
            elif strand_dict[gene_name] == -1:

                ### Get positions for gene start/end WITH UTR (-utr#, +utr#)
                pos = (location_dict[gene_name][0]-utr_length_to_include, location_dict[gene_name][1]+utr_length_to_include)
                if pos[1] < pos[0]:
                    print('found a bug: {} -- start-utr: {} end+utr: {} strand: {}'.format(gene_name, pos[0], pos[1], strand_dict[gene_name]))
                    continue

                ### Create a list of read values across the length of gene with UTRs; 0 if there are no reads at that position
                sequencing = []
                for i in range(pos[0]+1, pos[1]+1):
                    try:
                        sequencing.append(rev_dict[i])
                    except KeyError:
                        sequencing.append(0)

                feature_dict_meta[sample_name][gene_name] = sequencing[::-1] ### Note the CRUCIAL reverse here to get all reads in same direction

    return feature_dict_meta

def filter_genes(sample_names, feature_dict_meta, user_inputs, minCDS=18):

    print('Filter genes')
    print('length minimum: %snt CDS' % minCDS, 'coverage minimum: 20% nonzero reads in CDS', 'expression minimum: 0.5 average # reads in CDS')

    utr_length_to_include = user_inputs['utr_length_to_include']

    for sample in sample_names:
        
        print('##################################')
        print(sample, ' inititially has ', len(feature_dict_meta[sample]), 'genes (sanity check)')
        to_delete = [] 
        for gene in feature_dict_meta[sample]:
            
            ####Length cutoff is to make sure each CDS is >=100 nts long
            if len(feature_dict_meta[sample][gene]) < minCDS + (utr_length_to_include)*2:
                to_delete.append(gene)
                
            ####Coverage cutoff is to make sure that at least 20% of the values in the CDS are non-zero (specifically in coding region)
            elif np.percentile(feature_dict_meta[sample][gene][utr_length_to_include:-1*utr_length_to_include], 80) <= 0:
                to_delete.append(gene) ###Important parameter here, genes to discard based on coverage
    
            ####Finally cut the lowest expression genes (mean less than 0.5 in coding region)
            elif np.mean(feature_dict_meta[sample][gene][utr_length_to_include:-1*utr_length_to_include]) < 0.5:
                to_delete.append(gene)
    
        for gene in to_delete:
            feature_dict_meta[sample].pop(gene)
        print('now has', len(feature_dict_meta[sample]), 'after coverage and length threshold')

    return(feature_dict_meta)

def common_genes(sample_names, feature_dict_meta, genome_dict, user_inputs, RNAseq=False):
    
    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    
    
    # Create list/dicts of genes appearing in all datasets for a given condition
    genes_in_all_samples = []
    genes_per_condition = {}
    genes_per_conditionPlusWT = {}

    # list of genes in all samples
    for gene in genome_dict['full_sequence_dict'].keys():
        templist = []
        for sample in sample_names:
            if gene not in feature_dict_meta[sample]:
                templist.append(sample)
        if len(templist) == 0:
            genes_in_all_samples.append(gene)
    
    # list of common genes for a given experimental condition, specified in user_inputs
    for condition in user_inputs['condition_info']:
        genes_per_condition[condition] = []
        for gene in genome_dict['full_sequence_dict'].keys():
            templist = []
            for sample in user_inputs['condition_info'][condition]:
                if gene not in feature_dict_meta[sample]:
                    templist.append(sample)
            if len(templist) == 0:
                genes_per_condition[condition].append(gene)

    print('##################################\n')
    print('There are', len(genes_in_all_samples), 'genes that appear in all datasets')
    print('\n##################################\n')
    
    for condition in user_inputs['condition_info']:
        print('There are {} genes that appear in all {} samples'.format(len(genes_per_condition[condition]), condition))
    print('\n##################################\n')


    if RNAseq == True:

        expConditions = [condition for condition in user_inputs['condition_info'] if (condition not in ['WT', 'DMSO'] and 
                                                                                   user_inputs['condition_type'][condition] == 'RiboSeq')]
        for condition in expConditions:
            genes_per_conditionPlusWT[condition] = []
            for gene in genes_per_condition[condition]:
                if gene in genes_per_condition[WT] and (gene in genes_per_condition[user_inputs['rnaSeq_pairs'][WT]] and 
                                                        gene in genes_per_condition[user_inputs['rnaSeq_pairs'][condition]]):
                    genes_per_conditionPlusWT[condition].append(gene)

        for condition in expConditions:
            print('There are {} genes that appear in all {} + {} + {} + {} samples'.format(len(genes_per_conditionPlusWT[condition]), 
                                                                                           WT, user_inputs['rnaSeq_pairs'][WT], 
                                                                                           condition, user_inputs['rnaSeq_pairs'][condition]))
        print('\n##################################')
    
    else:
        expConditions = [condition  for condition in user_inputs['condition_info'] if condition not in ['WT', 'DMSO']]
        # list of common genes for a given experimental condition, combined with WT samples - this is useful for fold change analyses
        for condition in expConditions:
            genes_per_conditionPlusWT[condition] = []
            for gene in genes_per_condition[condition]:
                if gene in genes_per_condition[WT]:
                    genes_per_conditionPlusWT[condition].append(gene)


        for condition in expConditions:
            print('There are {} genes that appear in all {} + {} samples'.format(len(genes_per_conditionPlusWT[condition]), WT, condition))
        print('\n##################################')

    return genes_in_all_samples, genes_per_condition, genes_per_conditionPlusWT

def gene_lookup(geneName, genome_dict):
    for gene in genome_dict['full_sequence_dict'].keys():
        if gene.find(geneName) != -1:
            print(gene)
            return(gene)


def rewrite_wig(sample_files_dict, user_inputs, genome_dict, outdir):

    mapping_offset = user_inputs['mapping_offset']
    alignment_type = user_inputs['alignment_type']

    print('Rewrite Wig files to incorporate mapping offset %s' % mapping_offset)

    total_genome_len = genome_dict['total_genome_len']
    
    for sample_name in sample_files_dict:
        fwd, rev = sample_files_dict[sample_name]
        print('##### {}'.format(sample_name))
        
        ################################
        ### Get the adjusted fwd and rev genome coverage
        fwd_dict = {}
        rev_dict = {}
        
        

        ##### forward ######
        header = ''
        # creating a dictionary with all positions shifted by + mapping_offset
        with open(fwd) as infile:
            for line in enumerate(infile):
                # first TWO LINES = HEADER of the wig file
                if line[0] <= 1: 
                    header += line[1]
                else:
                    split_line = line[1].split('\t')   # wig format: [position, readCount]
                    adjusted_position = int(split_line[0]) + mapping_offset                                        # Note: mapping offset addition
                    readCount = float(split_line[1])
                    # wig is sorted by position, so read in low to high
                    # circular genome --> negative pos map to end of genome, pos > genome length map to start of genome
                    if adjusted_position < 0:                         # negative mapping offset only
                        fwd_dict[total_genome_len + adjusted_position] = readCount
                    elif (0 <= adjusted_position <= total_genome_len):
                        if adjusted_position in fwd_dict.keys():
                            fwd_dict[adjusted_position] += readCount
                        else:
                            fwd_dict[adjusted_position] = readCount
                    elif adjusted_position > total_genome_len:         # positive mapping offset only
                        if (adjusted_position - total_genome_len) in fwd_dict.keys():
                            fwd_dict[adjusted_position - total_genome_len] += readCount
                        else:
                            fwd_dict[adjusted_position - total_genome_len] = readCount
                            
        # sort dictionary by position
        fwd_writing = sorted(fwd_dict.items(), key=lambda x: x[0])
        # write to new wig file
        with open(outdir + sample_name + '_%s_fwd_offset_%s.wig' % (alignment_type, mapping_offset), 'w') as outfile:
            outfile.write(header)
            for i,j in fwd_writing:
                # adjusted_position (tab) readCount (enter)
                outfile.write('{}\t{}\n'.format(i,int(j)))
        print('Done with fwd')

        ##### reverse #####

        header = ''
        with open(rev) as infile:
            for line in enumerate(infile):
                # first TWO LINES = HEADER of the wig file
                if line[0] <= 1: 
                    header += line[1]
                else:     
                    split_line = line[1].split('\t')   # wig format: [position, readCount]
                    # positive mapping offset --> adj_pos can be negative or >= 0 and <total_len_genome
                    # negative mapping offset --> adj_pos can be > total_genome_len or >= 0 and <=total_len_genome
                    adjusted_position = int(split_line[0]) - mapping_offset                                        # Note: mapping offset subtraction
                    readCount = float(split_line[1])
                    # wig is sorted by position, so read in low to high
                    # circular genome --> negative pos map to end of genome, pos > genome length map to start of genome
                    if adjusted_position < 0:                         # positive mapping offset only
                        rev_dict[total_genome_len + adjusted_position] = readCount
                    elif (0 <= adjusted_position <= total_genome_len):
                        if adjusted_position in rev_dict.keys():
                            rev_dict[adjusted_position] += readCount
                        else:
                            rev_dict[adjusted_position] = readCount
                    elif adjusted_position > total_genome_len:         # negative mapping offset only
                        if (adjusted_position - total_genome_len) in rev_dict.keys():
                            rev_dict[adjusted_position - total_genome_len] += readCount
                        else:
                            rev_dict[adjusted_position - total_genome_len] = readCount

        # sort dictionary by position
        rev_writing = sorted(rev_dict.items(), key=lambda x: x[0])
        # write to new wig file
        with open(outdir + sample_name + '_%s_rev_offset_%s.wig' % (alignment_type, mapping_offset), 'w') as outfile:
            outfile.write(header)
            for i,j in rev_writing:
                # adjusted_position (tab) readCount (enter)
                outfile.write('{}\t{}\n'.format(i,int(j)))
        print('Done with rev')


def wig_to_cpm(sample_files_dict, total_read_dict, user_inputs, outdir):

    alignment_type = user_inputs['alignment_type']
    mapping_offset = user_inputs['mapping_offset']

    print('Wig --> _cpm.wig')
    
    #note: not adding mapping offset
    for sample_name in sample_files_dict:
        fwd, rev = sample_files_dict[sample_name]
        print('##### {}'.format(sample_name))

        ################################
        ### Get the adjusted fwd and rev genome coverage
        fwd_dict = {}
        rev_dict = {}

        
        ##### forward ######
        header = ''
        # creating a dictionary with all read counts converted to counts per million
        with open(fwd) as infile:
            for line in enumerate(infile):
                # first TWO LINES = HEADER of the wig file
                if line[0] <= 1: 
                    header += line[1]
                else:
                    split_line = line[1].split('\t')   # wig format: [position, readCount]
                    position = int(split_line[0])
                    readCount = float(split_line[1])
                    cpm = readCount * 10**6 / total_read_dict[sample_name]
                    fwd_dict[position] = cpm

        # sort dictionary by position
        fwd_writing = sorted(fwd_dict.items(), key=lambda x: x[0])
        # write to new wig file
        with open(outdir + sample_name + '_%s_fwd_offset_%s_cpm.wig' % (alignment_type, mapping_offset), 'w') as outfile:
            outfile.write(header)
            for i,j in fwd_writing:
                # position (tab) readCount (enter)
                outfile.write('{}\t{}\n'.format(i, j))
        print('Done with fwd')

        
        ##### reverse ######
        header = ''
        # creating a dictionary with all read counts converted to counts per million
        with open(rev) as infile:
            for line in enumerate(infile):
                # first TWO LINES = HEADER of the wig file
                if line[0] <= 1: 
                    header += line[1]
                else:
                    split_line = line[1].split('\t')   # wig format: [position, readCount]
                    position = int(split_line[0])
                    readCount = float(split_line[1])
                    cpm = readCount * 10**6 / total_read_dict[sample_name]
                    rev_dict[position] = cpm

        # sort dictionary by position
        rev_writing = sorted(rev_dict.items(), key=lambda x: x[0])
        # write to new wig file
        with open(outdir + sample_name + '_%s_rev_offset_%s_cpm.wig' % (alignment_type, mapping_offset), 'w') as outfile:
            outfile.write(header)
            for i,j in rev_writing:
                # position (tab) readCount (enter)
                outfile.write('{}\t{}\n'.format(i,j))
        print('Done with rev')

################ profiling analysis ###################################################################################################333


def geneRPKM(sample_names, feature_dict_meta, total_read_dict, genome_dict, user_inputs):
    
    print('gene-level RPKM calculation')
    #outputs a dataframe of gene x RPKM normalized readcount at the gene for a given sample
    df_master = pd.DataFrame()
    
    location_dict = genome_dict['location_dict']
    strand_dict = genome_dict['strand_dict']
    utr_length_to_include = user_inputs['utr_length_to_include']
    
    for sample in sample_names:
        print(sample)
        for gene in strand_dict.keys():
            
            # length = cds end - cds start
            geneLength = location_dict[gene][1] - location_dict[gene][0]
            
            # Reads across CDS only (no UTRs)
            reads = feature_dict_meta[sample][gene][utr_length_to_include : -1 * utr_length_to_include]
            
            # Calculate RPKM across CDS & place in datatable of gene vs sample
            # focusing on gene-level analysis only; codon analysis will be separate
            df_master.at[gene, sample+'_RPKM'] = ( np.sum(reads) * 10**9) / ( geneLength * total_read_dict[sample] )
    
            # Remove infinite values (will cause issues downstream) & replace with NaN
            df_master[sample + '_RPKM'].replace(np.inf, np.nan, inplace=True)
            df_master[sample+'_RPKM'].replace(-np.inf, np.nan, inplace=True)

    return(df_master)


def TE(df_RPKM, total_read_dict, user_inputs):
    
    print('Calculate TE for riboSeq-rnaSeq pairs:')
    for sample_pair in user_inputs['rnaSeq_pairs'].items():
        if sample_pair[0] in user_inputs['filenames']:
            print(sample_pair)
            ribo = sample_pair[0]
            rna = sample_pair[1]
            
            # create a new column for sample_TE & calculate RPKM ribo / RPKM rna
            df_RPKM[ribo + '_TE'] = df_RPKM[ribo + '_RPKM'] / df_RPKM[rna + '_RPKM']
            
            # Remove infinite values & convert to NaN
            df_RPKM[ribo + '_TE'].replace(np.inf, np.nan, inplace=True)
            df_RPKM[ribo + '_TE'].replace(-np.inf, np.nan, inplace=True)
        
    return(df_RPKM)

def conditionAvg(df_RPKM_TE, user_inputs, RNAseq = False):

    print('\nCalculate condition average')
    
    for condition in user_inputs['condition_info']:
        if user_inputs['condition_type'][condition] == 'RiboSeq':
            print(condition)
            for sample in user_inputs['condition_info'][condition]:
                print('\t',sample)
    
            if RNAseq == False:
                samples = [sample + '_RPKM' for sample in user_inputs['condition_info'][condition]]
                
                df_RPKM_TE[condition + '_avgRPKM'] = df_RPKM_TE[samples].mean(skipna=False, axis=1)
            else:
                samples = [sample + '_TE' for sample in user_inputs['condition_info'][condition] if user_inputs['condition_type'][condition] == 'RiboSeq']
            
                df_RPKM_TE[condition + '_avgTE'] = df_RPKM_TE[samples].mean(skipna=False, axis=1)
    
    return(df_RPKM_TE)



def logFoldChange(df_avg, user_inputs, RNAseq=False):
    
    print('\nCalculate log2 fold change compared to WT')

    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    
    for condition in user_inputs['condition_info']: 
        if user_inputs['condition_type'][condition] == 'RiboSeq':
            if condition != WT:
                if RNAseq == False:
                    df_avg[condition + '_RPKM_log2FC'] = df_avg[condition + '_avgRPKM'].apply(np.log2) - df_avg[WT + '_avgRPKM'].apply(np.log2)
                    #replace inf with NaN
                    df_avg[condition + '_RPKM_log2FC'].replace(np.inf, np.nan, inplace=True)
                    df_avg[condition + '_RPKM_log2FC'].replace(-np.inf, np.nan, inplace=True)
    
                else:
                    df_avg[condition + '_TE_log2FC'] = df_avg[condition + '_avgTE'].apply(np.log2) - df_avg[WT + '_avgTE'].apply(np.log2)
                    #replace inf with NaN
                    df_avg[condition + '_TE_log2FC'].replace(np.inf, np.nan, inplace=True)
                    df_avg[condition + '_TE_log2FC'].replace(-np.inf, np.nan, inplace=True)

    return(df_avg)


def independent_t_test(df_RPKM_TE, user_inputs, RNAseq):

    print('Assess statistical significance of differences in RPKM/TE using independent t-tests for each gene')
    print('reporting raw p-values')
    
    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    expConditions = [condition  for condition in user_inputs['condition_info'] if (condition != WT and user_inputs['condition_type'][condition] =='RiboSeq')]
    
    if RNAseq == False:
        datatype = 'RPKM'
    else:
        datatype = 'TE'

    
    # create an empty column to enter the pval info
    for condition in expConditions:
        df_RPKM_TE[condition + '_%s_rawTtest_pval' % datatype] = np.nan
        
    # Perform independent T-test for each gene for a given condition
    # compare all WT RPKM for that gene with all condition RPKM for that gene
    for gene in df_RPKM_TE.index: 
    
        a = []
        for sample in user_inputs['condition_info'][WT]:
            a.append(df_RPKM_TE.loc[gene][sample + '_%s' % datatype])

        for condition in expConditions:
            b = []
            for sample in user_inputs['condition_info'][condition]:
                b.append(df_RPKM_TE.loc[gene][sample + '_%s' % datatype])

            # assuming equal variance in the T-test -- may not be correct
            t_val, p_val = stats.ttest_ind(a, b, equal_var=True)
            df_RPKM_TE.at[gene, condition + '_%s_rawTtest_pval' % datatype] = p_val

    return df_RPKM_TE           
                 

def multipleTest(df_raw_pval, user_inputs, RNAseq=False):
    
    print('Multiple testing --> FDR correcting p-values')
    
    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    expConditions = [condition  for condition in user_inputs['condition_info'] if (condition != WT and user_inputs['condition_type'][condition] =='RiboSeq')]
    
    if RNAseq == False:
        datatype = 'RPKM'
    else:
        datatype = 'TE'
    
    for condition in expConditions:
        # create a temporary df where all genes with null p-values for that sample are removed
        temp_df = df_raw_pval[df_raw_pval[condition + '_%s_rawTtest_pval' % datatype].isnull() ==False]

        sigs, new_pvals, trash1, trash2 = statsmodels.stats.multitest.multipletests(list(temp_df[condition + '_%s_rawTtest_pval' % datatype]),\
                                                                                    alpha=0.05, method='fdr_bh')
        temp_df['fdr_corrected_pval_bh'] = new_pvals
        df_raw_pval[condition + '_%s_fdrBH_corrected_pval' % datatype] = np.nan
        for gene in temp_df.index:
            df_raw_pval.at[gene, condition + '_%s_fdrBH_corrected_pval' % datatype] = temp_df.loc[gene]['fdr_corrected_pval_bh']

    return df_raw_pval

        
########################################################## codon analysis ######################################################################

            
def codon_pauseScore(samples_to_consider, genes_to_consider, feature_dict_meta, genome_dict, total_read_dict, user_inputs, exclude_end_nt=30):

    print('calculate pause score by codon, excluding terminal %s nt and codons with < 5 reads' % exclude_end_nt)
    print('codon identifiers use base 1 positioning')
    print('########', ', '.join(samples_to_consider))

    condition_pauseScore_meta = {}
    condition_cpm_meta = {}
    condition_count_meta = {}

    location_dict = genome_dict['location_dict']
    utr_length_to_include = user_inputs['utr_length_to_include']
    
    for sample in samples_to_consider:
        condition_pauseScore_meta[sample] = {}
        condition_cpm_meta[sample] = {}
        condition_count_meta[sample] = {}
        
        for gene in genes_to_consider:
            
            # length = cds end - cds start
            geneLength = location_dict[gene][1] - location_dict[gene][0]

            # calculate reads across whole gene
            # CDS only (no UTRs)
            reads = feature_dict_meta[sample][gene][utr_length_to_include : -1 * utr_length_to_include]
            geneSum = np.sum(reads)

            # calculate reads across codon region - excluding 1st & last nt based on exclude_end_nt
            for i in range(exclude_end_nt, geneLength - exclude_end_nt):
                
                if i % 3 == 0 and i +2 <= geneLength:
                    codonSum = np.sum(reads[i: i+3])

                    # if codon has minimum of 5 reads, add to dictionary
                    if codonSum > 5:
                        pause_score = codonSum / geneSum *10**3
                        # codons named by geneName_codonNumber (base 1 for easier manual reference)
                        codon_id = gene + '_%s' % '{:04d}'.format(int(i/3) + 1)
    
                        condition_pauseScore_meta[sample][codon_id] = pause_score
                        condition_cpm_meta[sample][codon_id] = codonSum / total_read_dict[sample] * 10**6
                        condition_count_meta[sample][codon_id] = codonSum

    return condition_pauseScore_meta, condition_cpm_meta, condition_count_meta


def average_pauseScore(pauseScore_meta, codons_per_conditionPlusWT, user_inputs ):

    print('Calculate average pause score within condition')

    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    WTsamples = user_inputs['condition_info'][WT]
    
    avg_pauseScore_meta = {}
    
    for condition in pauseScore_meta:
        print('### ' + condition)
        avg_pauseScore_meta[condition] = {condition: {}, WT: {} }
        
        expSamples = user_inputs['condition_info'][condition]
        codons_to_consider = codons_per_conditionPlusWT[condition]
    
        for codon_id in codons_to_consider:
    
            # calculate average pause score across experimental condition
            exp_temp = [ pauseScore_meta[condition][sample][codon_id] for sample in expSamples]
            avg_pauseScore_meta[condition][condition][codon_id] = np.mean(exp_temp)
    
            # calculate average pause score across experimental condition
            wt_temp = [ pauseScore_meta[condition][sample][codon_id] for sample in WTsamples]
            avg_pauseScore_meta[condition][WT][codon_id] = np.mean(wt_temp)
            
    return avg_pauseScore_meta


def average_pauseScore_all(pauseScore_meta_all, codons_in_all_samples, user_inputs ):

    print('Calculate average pause score within condition, for all samples including WT/DMSO')
    
    avg_pauseScore_meta_all = {}
    
    for condition in user_inputs['condition_info']:
        print('### ' + condition)
        avg_pauseScore_meta_all[condition] = {}
        
        expSamples = user_inputs['condition_info'][condition]

        for codon_id in codons_in_all_samples:

            # calculate average pause score across experimental condition
            exp_temp = [ pauseScore_meta_all[sample][codon_id] for sample in expSamples]
            avg_pauseScore_meta_all[condition][codon_id] = np.mean(exp_temp)

    return avg_pauseScore_meta_all

    
        
def pauseScore_log2FC(avg_pauseScore_meta, user_inputs):

    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    
    print('Calculate fold change for condition vs %s, from average pause scores' % WT)

    log2FC_pauseScore_meta = {}
    
    for condition in avg_pauseScore_meta:
        print('### ' + condition)
        
        log2FC_pauseScore_meta[condition] = {}
        
        for codon_id in avg_pauseScore_meta[condition][WT]:
            log2FC = np.log2( avg_pauseScore_meta[condition][condition][codon_id] ) - np.log2( avg_pauseScore_meta[condition][WT][codon_id] ) 
    
            log2FC_pauseScore_meta[condition][codon_id] = log2FC
    
    return log2FC_pauseScore_meta
    

def pauseScore_log2FC_all(avg_pauseScore_meta_all, user_inputs):

    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    expConditions = [condition  for condition in user_inputs['condition_info'] if condition not in ['WT', 'DMSO']]
    
    print('Calculate fold change for condition vs %s, from average pause scores - all samples' % WT)

    log2FC_pauseScore_meta_all = {}

    for condition in expConditions:
        print('### ' + condition)

        log2FC_pauseScore_meta_all[condition] = {}
        
        for codon_id in avg_pauseScore_meta_all[WT]:
            log2FC = np.log2( avg_pauseScore_meta_all[condition][codon_id] ) - np.log2( avg_pauseScore_meta_all[WT][codon_id] )
            
            log2FC_pauseScore_meta_all[condition][codon_id] = log2FC

    return log2FC_pauseScore_meta_all
    
   

def get_peak_seq(codon_dict, outfile, genome_dict, user_inputs, seqWindow=(27,6)):
    
    full_sequence_dict = genome_dict['full_sequence_dict']
    location_dict = genome_dict['location_dict']
    utr_length_to_include = user_inputs['utr_length_to_include']

    aaSeq_list = []

    for codon_id in codon_dict.keys():
        
        #codon_id_meta[sample][gene][int(i/3)] = gene+'_%s' % '{:04d}'.format(int(i/3) + 1)
        codon_id = codon_id.split('_')
        gene_name = '_'.join(codon_id[:-1])
        codon_num = int(codon_id[-1])
        
        #reverse the base 1 positioning to base 0
        nt_pos = (codon_num - 1) *3
        start, stop = location_dict[gene_name]
        
        geneSeq = full_sequence_dict[gene_name][utr_length_to_include:-utr_length_to_include]
        # seqWindow (27, 6) gives 11aa seq ending at A site, for mapping_offset aligned to P site
        if nt_pos - seqWindow[0]>= 0:
            peak_seq = geneSeq[nt_pos-seqWindow[0]:nt_pos + seqWindow[1]]
        else: 
            peak_seq = geneSeq[0:nt_pos + seqWindow[1]]
        #print(codon_id, geneSeq, peak_seq)
        
        aaSeq_list.append(str(Seq.translate(peak_seq)))
        
    with open(outfile, 'w') as f:
        f.write('\n'.join(aaSeq_list))

    return aaSeq_list



def write_pauseScoreCSV(sorted_log2FC_pauseScore, aaSeq_dict, avg_pauseScore_meta, codon_cpm_meta, pauseScore_meta, genome_dict, outdir, user_inputs):

    location_dict = genome_dict['location_dict']
    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    
    for condition in sorted_log2FC_pauseScore:
        samplelist = list(codon_cpm_meta[condition].keys())
        conditionlist = list(avg_pauseScore_meta[condition].keys())
        print(samplelist, conditionlist)
        codon_stalls = [ ]
        
        for entry in map(list, zip(sorted_log2FC_pauseScore[condition].keys(),  aaSeq_dict[condition]['all'], sorted_log2FC_pauseScore[condition].values())):
            codon_stalls.append(entry)
            
        for i in range(0, len(codon_stalls)):
            codon_id = codon_stalls[i][0]
            codon_split = codon_id.split('_')
            gene_id = '_'.join(codon_split[:-1])
            codon_num = int(codon_split[-1])
            geneLength = location_dict[gene_id][1] - location_dict[gene_id][0]

            codon_stalls[i].insert(1, gene_id)
            codon_stalls[i].insert(2, codon_num)
            
            # append codon cpm for all samples in samplelist
            for sample in samplelist:
                codon_stalls[i] += [codon_cpm_meta[condition][sample][codon_id]]
                
            # append codon rpkm for all samples in samplelist
            for sample in samplelist:
                codon_stalls[i] += [codon_cpm_meta[condition][sample][codon_id] *10**3 / geneLength]

            # append codon pause score for all samples in samplelist
            for sample in samplelist:
                codon_stalls[i] += [pauseScore_meta[condition][sample][codon_id]]

            # append codon average pause score for each condition
            for sample in conditionlist:
                codon_stalls[i] += [avg_pauseScore_meta[condition][sample][codon_id]]

    
        header = ['codonID', 'geneID', 'codonNum (base1)', 'aaSeq', 'log2FC %s/%s' % (condition, WT)]
        header += [x + '_cpm' for x in samplelist] + [x + '_rpkm' for x in samplelist] 
        header += [x + '_pauseScore' for x in samplelist] + [x + '_avg_pauseScore' for x in conditionlist] 
        print(header, '\n', codon_stalls[0], '\n')
        
        
        # export to csv
        with open(outdir + condition + '_codonStalls.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(header)
            writer.writerows(codon_stalls)




def get_codon_aa_freq(sorted_log2FC_pauseScore_all, genome_dict, user_inputs):

    full_sequence_dict = genome_dict['full_sequence_dict']
    location_dict = genome_dict['location_dict']
    utr_length_to_include = user_inputs['utr_length_to_include']
    
    codon_freq = { 'E': {}, 'P': {}, 'A': {} }
    aa_freq = { 'E': {}, 'P': {}, 'A': {} }

    for condition in sorted_log2FC_pauseScore_all:
            
        codon_freq['E'][condition] = {}
        codon_freq['P'][condition] = {}
        codon_freq['A'][condition] = {}

        aa_freq['E'][condition] = {}
        aa_freq['P'][condition] = {}
        aa_freq['A'][condition] = {}

        for codon in gencode.keys():
        
            codon_freq['E'][condition][codon] = []
            codon_freq['P'][condition][codon] = []
            codon_freq['A'][condition][codon] = []

        for aa in aaList:
        
            aa_freq['E'][condition][aa] = []
            aa_freq['P'][condition][aa] = []
            aa_freq['A'][condition][aa] = []

        
    for condition in sorted_log2FC_pauseScore_all:

        for codon_id in sorted_log2FC_pauseScore_all[condition]:
            
            #codon_id  = gene+'_%s' % '{:04d}'.format(int(i/3) + 1)
            codon_det = codon_id.split('_')
            gene_name = '_'.join(codon_det[:-1])
            codon_num = int(codon_det[-1])

            #reverse the base 1 positioning to base 0
            nt_pos = (codon_num - 1) *3
            start, stop = location_dict[gene_name]

            geneSeq = full_sequence_dict[gene_name][utr_length_to_include:-utr_length_to_include]
            
            E = geneSeq[nt_pos-3:nt_pos]
            P = geneSeq[nt_pos:nt_pos+3]
            A = geneSeq[nt_pos+3:nt_pos+6]

            codon_freq['E'][condition][E].append(sorted_log2FC_pauseScore_all[condition][codon_id])
            codon_freq['P'][condition][P].append(sorted_log2FC_pauseScore_all[condition][codon_id])
            codon_freq['A'][condition][A].append(sorted_log2FC_pauseScore_all[condition][codon_id])

            aa_freq['E'][condition][Seq.translate(E)].append(sorted_log2FC_pauseScore_all[condition][codon_id])
            aa_freq['P'][condition][Seq.translate(P)].append(sorted_log2FC_pauseScore_all[condition][codon_id])
            aa_freq['A'][condition][Seq.translate(A)].append(sorted_log2FC_pauseScore_all[condition][codon_id])

        for codon in gencode.keys():
            
            avg_log2FC_codon = np.mean(codon_freq['E'][condition][codon])
            codon_freq['E'][condition][codon] = avg_log2FC_codon

            avg_log2FC_codon = np.mean(codon_freq['P'][condition][codon])
            codon_freq['P'][condition][codon] = avg_log2FC_codon

            avg_log2FC_codon = np.mean(codon_freq['A'][condition][codon])
            codon_freq['A'][condition][codon] = avg_log2FC_codon

        for aa in aaList:

            avg_log2FC_aa = np.mean(aa_freq['E'][condition][aa])
            aa_freq['E'][condition][aa] = avg_log2FC_aa

            avg_log2FC_aa = np.mean(aa_freq['P'][condition][aa])
            aa_freq['P'][condition][aa] = avg_log2FC_aa

            avg_log2FC_aa = np.mean(aa_freq['A'][condition][aa])
            aa_freq['A'][condition][aa] = avg_log2FC_aa
    
    return codon_freq, aa_freq



##########################################################################################################################
def get_peak_seq_single(codon_id, genome_dict, user_inputs, seqWindow=(27,6)):
    
    full_sequence_dict = genome_dict['full_sequence_dict']
    location_dict = genome_dict['location_dict']
    utr_length_to_include = user_inputs['utr_length_to_include']

    #codon_id_meta[sample][gene][int(i/3)] = gene+'_%s' % '{:04d}'.format(int(i/3) + 1)
    codon_id = codon_id.split('_')
    gene_name = '_'.join(codon_id[:-1])
    codon_num = int(codon_id[-1])
    
    #reverse the base 1 positioning to base 0
    nt_pos = (codon_num - 1) *3
    start, stop = location_dict[gene_name]
    
    geneSeq = full_sequence_dict[gene_name][utr_length_to_include:-utr_length_to_include]
    # seqWindow (27, 6) gives 11aa seq ending at A site, for mapping_offset aligned to P site
    if nt_pos - seqWindow[0]>= 0:
        peak_seq = geneSeq[nt_pos-seqWindow[0]:nt_pos + seqWindow[1]]
    else: 
        peak_seq = geneSeq[0:nt_pos + seqWindow[1]]
    #print(codon_id, geneSeq, peak_seq)
        
    return str(Seq.translate(peak_seq))
