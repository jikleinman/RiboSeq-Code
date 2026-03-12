########## imports #############

from Bio import SeqIO, Seq
import numpy as np
from collections import Counter
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import cm, colors
from matplotlib.ticker import MultipleLocator
import copy
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
import statsmodels.stats
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import RiboSeq_Analysis_fxns as analysis


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



###########################################################################################################################################
########################################## PLOT SETTINGS ##################################################################################
###########################################################################################################################################


matplotlib.rcParams['xtick.labelsize'] = 10
matplotlib.rcParams['ytick.labelsize'] = 10
matplotlib.rcParams['axes.labelsize'] = 10
matplotlib.rcParams['axes.titlesize'] = 10

matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['grid.color'] = '0.5'
matplotlib.rcParams['grid.linewidth'] = '0.5'

matplotlib.rcParams['axes.edgecolor'] = '0.25'
matplotlib.rcParams['xtick.color'] = '0'
matplotlib.rcParams['ytick.color'] = '0'

matplotlib.rcParams['xtick.major.width'] = 1
matplotlib.rcParams['ytick.major.width'] = 1
matplotlib.rcParams['ytick.major.size'] = 5
matplotlib.rcParams['xtick.major.size'] = 5
matplotlib.rcParams['axes.spines.right'] = True
matplotlib.rcParams['axes.spines.left'] = True
matplotlib.rcParams['axes.spines.top'] = True
matplotlib.rcParams['axes.spines.bottom'] = True
matplotlib.rcParams['axes.axisbelow'] = True
matplotlib.rcParams['figure.dpi'] = 600

# add colorbar to bottom for stall strength
BuViPi = LinearSegmentedColormap.from_list('BuViPi', (
                    # Edit this gradient at https://eltos.github.io/gradient/#DC267E-795EF0-56B4E9
                    (0.000, (0.863, 0.149, 0.494)),
                    (0.500, (0.475, 0.369, 0.941)),
                    (1.000, (0.337, 0.706, 0.914))))

BuWhPi = LinearSegmentedColormap.from_list('BuWhPi', (
                    # Edit this gradient at https://eltos.github.io/gradient/#DC267E-795EF0-56B4E9
                    (0.000, (0.863, 0.149, 0.494)),
                    (0.500, (1.000, 1.000, 1.000)),
                    (1.000, (0.337, 0.706, 0.914))))

ViWhPi = LinearSegmentedColormap.from_list('ViWhiPi', (
    # Edit this gradient at https://eltos.github.io/gradient/#ViWhiPi=DC267E-FFFFFF-795EF0
    (0.000, (0.863, 0.149, 0.494)),
    (0.500, (1.000, 1.000, 1.000)),
    (1.000, (0.475, 0.369, 0.941))))

###########################################################################################################################################
########################################## METAGENE PLOTS #################################################################################
###########################################################################################################################################


def plot_periodicity(sample_names, feature_dict_meta, outdir, user_inputs):
    
    mapping_offset = user_inputs['mapping_offset']
    utr_length_to_include = user_inputs['utr_length_to_include']
    print('Plot periodicity in CDS regions')
    print('excluding 1st and last 18nt of CDS & filtering to remove low expression/coverage genes (mean # reads <  0.5 or > 80% reads = 0)')
    for sample in sample_names:
        print('####### {}'.format(sample))

        # count number of genes passing filters (avg reads/position in gene > 0.5, < 80% reads/position = 0 ) 
        valid_genes = 0
        # create a list for counting normalized reads at each position in codon across whole genome
        first, second, third = [], [], []
        for name, reads in feature_dict_meta[sample].items():
            
            ### Extract just the CDS
            cds_reads = reads[utr_length_to_include:-1*utr_length_to_include]
            ### Take only the center of the CDS to remove initiation/termination effects
            cds_reads = cds_reads[18:-18]

            # get the average read count across the gene, with start/stop excluded
            meany = np.mean(cds_reads)
            
            ### Filter out & exclude low expression/coverage genes from this analysis
            if meany <= 0.5: 
                continue                                   #exclude genes with average <0.5 reads/position
            if np.percentile(cds_reads, 80) <= 0:
                continue                                   # exclude genes where >= 80% reads are 0 (?)
            '''if len(cds_reads) < 100 -18*2:                 # exclude genes where length after start/stop removal is < 64 nt (21 aa)
               continue'''
                
            ### Normalize reads for each gene -- read count by position / avg read counts for whole gene (start/stop excluded)
            cds_reads = np.array(cds_reads)/meany
            ### add normalized read counts for each position in gene to a list, separated by position in codon 
            temp_first = cds_reads[0::3] # starting from position zero & taking every 3rd val
            temp_second = cds_reads[1::3]
            temp_third = cds_reads[2::3]
            first.extend(temp_first)
            second.extend(temp_second)
            third.extend(temp_third)
            valid_genes += 1

        # plot periodicity via total normalized read counts per codon position across all CDS regions specified
        fig, ax = plt.subplots()
        ax.bar([1, 2, 3], [np.mean(first), np.mean(second), np.mean(third)],\
                   yerr=[np.std(first)/len(first), np.std(second)/len(second), np.std(third)/len(third)])
        ax.set_title('{}, n={} genes'.format(sample, valid_genes))
        ax.set_ylabel('Average relative abundance')
        ax.set_xlabel('Codon positions')
        ax.set_ylim(0, 1.5)
        ax.set_xticks([1,2,3])
        #plt.show()
        plt.savefig(outdir + sample + '_read_periodicity_offset_%s.pdf' % mapping_offset, bbox_inches='tight')




### 3' and 5' metagene analyses are normalized by average reads across the gene & total sequencing depth 
# read / (avg read value for gene * total sequencing depth) * 10^9
import math
def plot_fiveprime_meta(feature_dict_meta, samples_to_consider, genes_to_consider, outdir, user_inputs, total_read_dict, analysis_type='total'):

    print(samples_to_consider)
    mapping_offset = user_inputs['mapping_offset']
    utr_length_to_include = user_inputs['utr_length_to_include']

    print('5 prime Metagene analysis')
    print(list(samples_to_consider), '%s genes,' % len(genes_to_consider), 'offset %s' % mapping_offset)
    print('normalize to average reads by gene & total sequencing depth')
    
    plotting_dict = {}
    full_dict_temp = {}
    winsorize = True

    ### get read counts around the start codon (1/2 UTR before start & 3 UTR-length after) --> normalize to mean reads for gene & total sequencing depth
    for sample in samples_to_consider:
        temp_array = []
        for gene in genes_to_consider:
            if winsorize:
                reads = stats.mstats.winsorize(np.array(feature_dict_meta[sample][gene]), axis=0, limits=0.05)
            else:
                reads = feature_dict_meta[sample][gene]
            # Calculate the mean reads across the gene, based on CDS only
            meany = np.mean(reads[utr_length_to_include:-1*utr_length_to_include]) 
             # normalize reads/position to mean reads for gene in coding region & to total sequencing depth of the sample
            temp_array.append([i/meany/int(total_read_dict[sample])*10**9 for i in reads[int(utr_length_to_include/2):utr_length_to_include*3]])
            
        temp_array = stats.mstats.winsorize(np.array(temp_array), axis=0, limits=0.05)
        
        y_vals = np.mean(np.array(temp_array), axis=0)
        plotting_dict[sample] = y_vals
        full_dict_temp[sample] = temp_array
       
    
    
    # list of x axis values
    x_vals = np.arange(-int(utr_length_to_include/2),utr_length_to_include*2, 1)
    
    ### plot normalized reads surrounding the start codon of all communal genes
    fig, ax = plt.subplots(figsize=(10,6))
    for sample in samples_to_consider:
        ax.plot(x_vals, plotting_dict[sample], alpha=0.8, label=sample)
    
    lgd = ax.legend(fontsize=16, loc='center left', bbox_to_anchor=(1, 0.5))
    
    
    ax.set_title('{} genes appearing in all {} {} datasets'.format(len(temp_array), len(samples_to_consider), analysis_type), fontsize=20)
    ax.tick_params(labelsize=16)
    ax.set_xlabel('Position relative to start', fontsize=20)
    ax.set_ylabel('Avg normalized reads', fontsize=20)
    maxy = max([max(plotting_dict[i]) for i in samples_to_consider]) 
    plt.ylim(top=175)
    ax.xaxis.set_major_locator(MultipleLocator(25))
    ax.axvline(0, c='k')
    plt.text(0,maxy+.1, '-1st nt of start codon')
    ax.axvline(3, c='k', linestyle = ':')
    plt.text(3,1, '+3')
    #ax.axvline(1,c ='r', linestyle =':')
    #ax.axvline(2,c ='r', linestyle =':')
    plt.tight_layout()
    if winsorize:
        plt.savefig('{}/5prime_{}_offset_{}_winsor.pdf'.format(outdir, analysis_type, mapping_offset), bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        plt.show()
        plt.savefig('{}/5prime_{}_offset_{}_nowinsor.pdf'.format(outdir, analysis_type, mapping_offset), bbox_inches='tight')

    
    ### calculate p-values over the same window
    print('calculating p-value from mixed linear model')

    pval_listy = []
    for position in range(len(x_vals)):
        df_stat = pd.DataFrame()
        
        for sample, temp_array in full_dict_temp.items():
            
            sample_name = sample
            # temp_array = [[normalized reads for gene1], [normalized reads for gene2],...]
            # vals = transposed temp_array at a given position along gene = [normalized read for gene1 at position0, normalized read for gene2 at position0,...]
            vals = np.array(temp_array).T[position]

            #create a dataframe with 3 lists & where values = list of normalized reads for all genes at a single position
            temp_df = pd.DataFrame()
            temp_df['value'] = vals
            temp_df['Condition'] = sample_name
            temp_df['gene_n'] = list(range(len(vals)))
            
            # append temp_df to the larger stats dataframe -- all samples in 1 df for a single position
            df_stat = pd.concat([df_stat, temp_df])
            
        # calculate p-values for the single position, grouping reads from all samples for each individual gene 
        # Mixed models allow accounting for of multiple sources of variability simultaneously - I think this is sample variation & differences by gene
        md = smf.mixedlm('value~C(Condition)', df_stat, groups=df_stat['gene_n']).fit()
        # add pvalues to list
        pval_listy.append(md.pvalues.iloc[1])
        
        
        

    ### plot p-values (x = position along gene, y = log(mixed linear model p-values) )
    fig, ax = plt.subplots(figsize=(8,6))
    # convert p-values to log scale
    ax.semilogy(x_vals, pval_listy)
    ax.axhline(0.01/100, c='r', linestyle='--')
    ax.tick_params(labelsize=16)
    ax.set_xlabel('Position relative to start', fontsize=20)
    ax.set_ylabel('p-value from mixedlm', fontsize=20)
    ax.axvline(0, c='k')
    ax.axvline(3, c='k', linestyle = ':')
    plt.text(3,1, '+3')
    ax.xaxis.set_major_locator(MultipleLocator(25))
    plt.tight_layout()
    if winsorize:
        plt.savefig('{}/5prime_{}_offset{}_pval_winsor.pdf'.format(outdir, analysis_type, mapping_offset), bbox_inches='tight')
    else:
        plt.savefig('{}/5prime_{}_offset{}_pval_nowinsor.pdf'.format(outdir, analysis_type, mapping_offset), bbox_inches='tight')
        



def plot_threeprime_meta(feature_dict_meta, samples_to_consider, genes_to_consider, outdir, user_inputs, total_read_dict, analysis_type='total'):

    mapping_offset = user_inputs['mapping_offset']
    utr_length_to_include = user_inputs['utr_length_to_include']

    print('3 prime Metagene analysis')
    print(list(samples_to_consider), '%s genes,' % len(genes_to_consider), 'offset %s' % mapping_offset)
    print('normalize to average reads by gene & total sequencing depth')
    
    plotting_dict = {}
    full_dict_temp = {}
    winsorize = True


    

    ### get read counts around the start codon (1/2 UTR before start & 3 UTR-length after) --> normalize to mean reads for gene & total sequencing depth
    for sample in samples_to_consider:
        temp_array = []
        for gene in genes_to_consider:
            if winsorize:
                reads = stats.mstats.winsorize(np.array(feature_dict_meta[sample][gene]), axis=0, limits=0.05)
            else:
                reads = feature_dict_meta[sample][gene]
            # Calculate the mean reads across the gene, based on CDS only
            meany = np.mean(reads[utr_length_to_include:-1*utr_length_to_include]) 
             # normalize reads/position to mean reads for gene in coding region & to total sequencing depth of the sample
            temp_array.append([i/meany/int(total_read_dict[sample])*10**9 for i in reads[-1*utr_length_to_include*3:-int(utr_length_to_include/2)]])
            
            
        temp_array = stats.mstats.winsorize(np.array(temp_array), axis=0, limits=0.05)
        
        y_vals = np.mean(np.array(temp_array), axis=0)
        plotting_dict[sample] = y_vals
        full_dict_temp[sample] = temp_array
       
    
    
    # list of x axis values
    x_vals = np.arange(-1*utr_length_to_include*2, int(utr_length_to_include/2), 1)
    
    ### plot normalized reads surrounding the stop codon of all communal genes    
    fig, ax = plt.subplots(figsize=(10,6))
    for sample in samples_to_consider:
        ax.plot(x_vals, plotting_dict[sample], alpha=0.8, label=sample)
    
    lgd = ax.legend(fontsize=16, loc='center left', bbox_to_anchor=(1, 0.5))
    
    
    ax.set_title('{} genes appearing in all {} {} datasets'.format(len(temp_array), len(samples_to_consider), analysis_type), fontsize=20)
    ax.tick_params(labelsize=16)
    ax.set_xlabel('Position relative to stop', fontsize=20)
    ax.set_ylabel('Avg normalized reads', fontsize=20)
    maxy = max([max(plotting_dict[i]) for i in samples_to_consider])
    plt.ylim(top=175)
    ax.xaxis.set_major_locator(MultipleLocator(25))
    ax.axvline(-3, c='k')
    plt.text(-3,maxy+.1, '-1st nt of stop codon (-3)')
    ax.axvline(-6, c='k', linestyle = ':')
    plt.text(-9, maxy+.1, '-6')
    #ax.axvline(-2,c ='r', linestyle =':')
    #ax.axvline(-1,c ='r', linestyle =':')
    
    plt.tight_layout()
    if winsorize:
        plt.savefig('{}/3prime_{}_offset_{}_winsor.pdf'.format(outdir, analysis_type, mapping_offset), bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        plt.savefig('{}/3prime_{}_offset_{}_nowinsor.pdf'.format(outdir, analysis_type, mapping_offset), bbox_inches='tight')
        

    ### calculate p-values over the same window
    print('calculating p-value from mixed linear model')

    pval_listy = []
    for position in range(len(x_vals)):
        df_stat = pd.DataFrame()
        
        for sample, temp_array in full_dict_temp.items():
            
            sample_name = sample
            # temp_array = [[normalized reads for gene1], [normalized reads for gene2],...]
            # vals = transposed temp_array at a given position along gene = [normalized read for gene1 at position0, normalized read for gene2 at position0,...]
            vals = np.array(temp_array).T[position]

            #create a dataframe with 3 lists & where values = list of normalized reads for all genes at a single position
            temp_df = pd.DataFrame()
            temp_df['value'] = vals
            temp_df['Condition'] = sample_name
            temp_df['gene_n'] = list(range(len(vals)))
            
            # append temp_df to the larger stats dataframe -- all samples in 1 df for a single position
            df_stat = pd.concat([df_stat, temp_df])
            
        # calculate p-values for the single position, grouping reads from all samples for each individual gene 
        # Mixed models allow accounting for of multiple sources of variability simultaneously - I think this is sample variation & differences by gene
        md = smf.mixedlm('value~C(Condition)', df_stat, groups=df_stat['gene_n']).fit()
        # add pvalues to list
        pval_listy.append(md.pvalues.iloc[1])
        #print(x_vals[position], md.pvalues.iloc[1], (math.log10(md.pvalues.iloc[1])))

    ### plot p-values (x = position along gene, y = log(mixed linear model p-values) )
    fig, ax = plt.subplots(figsize=(8,6))
    # convert p-values to log scale
    ax.semilogy(x_vals, pval_listy)
    ax.axhline(0.01/100, c='r', linestyle='--')
    ax.tick_params(labelsize=16)
    ax.set_xlabel('Position relative to stop', fontsize=20)
    ax.set_ylabel('p-value from mixedlm', fontsize=20)
    miny = min(pval_listy)
    ax.axvline(-3, c='k')
    plt.text(-3, miny *5, '-3')
    ax.axvline(-6, c='k', linestyle = ':')
    plt.text(-9 ,miny *5, '-6')
    ax.xaxis.set_major_locator(MultipleLocator(25))
    plt.tight_layout()
    if winsorize:
        plt.savefig('{}/3prime_{}_offset{}_pval_winsor.pdf'.format(outdir, analysis_type, mapping_offset), bbox_inches='tight')
    else:
        plt.savefig('{}/3prime_{}_offset{}_pval_nowinsor.pdf'.format(outdir, analysis_type, mapping_offset), bbox_inches='tight')



# positional metagene = redistribution of reads across full gene window

def lengthwise_metagene(samples_to_consider, genes_to_consider, feature_dict_meta, outdir, user_inputs, analysis_type='total'):  
    mapping_offset = user_inputs['mapping_offset']
    utr_length_to_include = user_inputs['utr_length_to_include']
    
    # positions by fraction of gene length, to 2nd decimal place. Base 0.01
    vals = [i/100 for i in range(1,101)]
    fig, ax = plt.subplots()

    if analysis_type == 'RiboSeq':
        samples = [s for s in samples_to_consider if s in user_inputs['rnaSeq_pairs'].keys()]
    elif analysis_type == 'RNAseq':
        samples = [s for s in samples_to_consider if s in user_inputs['rnaSeq_pairs'].values()]
    else:
        samples = samples_to_consider
    
    for sample in samples:
        
        xvals = {i:0 for i in vals}
        for gene in genes_to_consider:
            reads = feature_dict_meta[sample][gene][utr_length_to_include:-utr_length_to_include]
            T = sum(reads)
            for i in range(0, len(reads)):
                pos = round((i+1)/len(reads), 2)
                if pos == 0:
                    pos=0.01
                # reads at a given position divided by total number of reads for that gene
                xvals[pos] += reads[i]/T
        for pos in xvals:
            # normalize to number of genes being considered to get average fraction of total gene reads
            xvals[pos] = xvals[pos]/len(genes_to_consider)
        plt.plot(xvals.keys(), xvals.values(), label=sample)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    ax.set_title('Lengthwise Metagene')
    ax.set_xlabel('Fraction of gene length')
    ax.set_ylabel('Reads per Total Gene Reads')
    plt.tight_layout()
    plt.savefig('{}/lengthwiseMetagene_{}_offset{}.pdf'.format(outdir, analysis_type, mapping_offset), bbox_inches='tight')

    


def short_CDS_window_rpm(geneID, samples_to_consider, total_read_dict, feature_dict_meta, genome_dict, outdir, user_inputs, window =[0,36]):
    if window[0] <0:
        window[0] = 0
    print('Plot gene window %s-%s - RPM normalized - %s' % (window[0], window[1], geneID))
    
    utr_length_to_include = user_inputs['utr_length_to_include']
    location_dict = genome_dict['location_dict']
    plotting_dict = {}
    
    fig,ax = plt.subplots()
    
    for sample in samples_to_consider:
        
        temp_array = []
        
        reads = feature_dict_meta[sample][geneID][utr_length_to_include:-utr_length_to_include]
        reads = reads[window[0]: window[1]]

        # calculate RPM at each location
        temp_array.append([i/int(total_read_dict[sample])*10**6 for i in reads])

        # y axis values    
        # since temp_array is a single list, np.mean on axis 0 gives an array of the same values but with appropriate dimensions
        y_vals = np.mean(np.array(temp_array), axis=0)
        plotting_dict[sample] = y_vals

    # list of x axis values, same length as reads
    x_vals = np.arange(window[0], window[1], 1)
    
    for sample in samples_to_consider:
        
        ax.plot(x_vals, plotting_dict[sample], alpha=0.8, label=sample)

    
    seq = genome_dict['full_sequence_dict'][geneID][utr_length_to_include:-utr_length_to_include]
    aaSeq = str(Seq.translate(seq))

    # label with aa & tri-nt sequence at every 3rd nt
    x_ticks = [x for x in range(window[0], window[1]) if x % 3==0]
    x_labels = []
    for i in range(window[0], window[1]):
        if i % 3 == 0:
            x_labels.append(str(int(i/3)+1) +  '\n' + aaSeq[int(i/3)] + '\n' + str(seq[i:i+3])) 

    ax.legend(fontsize=16, loc='center left', bbox_to_anchor=(1, 0.5))
    
    ax.set_ylabel('reads per million (RPM)', fontsize=16)    
    ax.set_xlabel(geneID.split('_')[2], style = 'italic', fontsize=16)
    ax.set_title('{} (aa {}-{})'.format(geneID.split('_')[2], int(window[0]/3)+1, int(window[1]/3)), fontsize=20)
    ax.set_xticks(x_ticks, x_labels)
    plt.savefig('{}/{}_rpm_aa-{}-{}.pdf'.format(outdir, geneID.split('_')[2], int(window[0]/3)+1, int(window[1]/3)), bbox_inches='tight')



def short_CDS_window_rpkm(geneID, samples_to_consider, total_read_dict, feature_dict_meta, genome_dict, outdir, user_inputs, window =[0,36]):
    
    print('Plot gene window %s-%s - RPKM normalized - %s' % (window[0], window[1], geneID))
    
    utr_length_to_include = user_inputs['utr_length_to_include']
    location_dict = genome_dict['location_dict']
    plotting_dict = {}
    
    fig,ax = plt.subplots()
    
    for sample in samples_to_consider:
        
        temp_array = []
        
        cds_reads = feature_dict_meta[sample][geneID][utr_length_to_include:-utr_length_to_include]
        reads = cds_reads[window[0]: window[1]]
        geneLength = location_dict[geneID][1] - location_dict[geneID][0]
        
        # calculate rpkm at each position
        temp_array.append([i/geneLength/int(total_read_dict[sample])*10**9 for i in reads])

        # y axis values    
        # since temp_array is a single list, np.mean on axis 0 gives an array of the same values but with appropriate dimensions
        y_vals = np.mean(np.array(temp_array), axis=0)
        plotting_dict[sample] = y_vals

    # list of x axis values, same length as reads
    x_vals = np.arange(window[0], window[1], 1)
    
    for sample in samples_to_consider:
        
        ax.plot(x_vals, plotting_dict[sample], alpha=0.8, label=sample)

    
    seq = genome_dict['full_sequence_dict'][geneID][utr_length_to_include:-utr_length_to_include]
    aaSeq = str(Seq.translate(seq))

    # label with aa & tri-nt sequence at every 3rd nt
    x_ticks = [x for x in range(window[0], window[1]) if x % 3==0]
    x_labels = []
    for i in range(window[0], window[1]):
        if i % 3 == 0:
            x_labels.append(str(int(i/3)+1) +  '\n' + aaSeq[int(i/3)] + '\n' + str(seq[i:i+3])) 
    
    ax.legend(fontsize=16, loc='center left', bbox_to_anchor=(1, 0.5))
    
    ax.set_ylabel('reads per kilobase per million (RPKM)', fontsize=16)    
    ax.set_xlabel(geneID.split('_')[2], style = 'italic', fontsize=16)
    ax.set_title('{} (aa {}-{})'.format(geneID.split('_')[2], int(window[0]/3)+1, int(window[1]/3)), fontsize=20)
    ax.set_xticks(x_ticks, x_labels)
    plt.savefig('{}/{}_rpkm_aa-{}-{}.pdf'.format(outdir, geneID.split('_')[2], int(window[0]/3)+1, int(window[1]/3)), bbox_inches='tight')


def short_CDS_window_TE(geneID, samples_to_consider, total_read_dict, feature_dict_meta, genome_dict, outdir, user_inputs, window =[0,36], RNAseq=False):
    
    if RNAseq == False:
        return(short_CDS_window_rpm(geneID, samples_to_consider, total_read_dict, feature_dict_meta, genome_dict, outdir, user_inputs, window))
    else:
        print('Plot gene window %s-%s - TE normalized (RPKM Ribo/ RPKM RNA) - %s' % (window[0], window[1], geneID))
    
    utr_length_to_include = user_inputs['utr_length_to_include']
    location_dict = genome_dict['location_dict']
    plotting_dict = {}
    
    fig,ax = plt.subplots()
    
    for sample in samples_to_consider:
        if sample in user_inputs['rnaSeq_pairs'].keys():
            rn = user_inputs['rnaSeq_pairs'][sample]
            
            temp_array = []
            
            sample_cds_reads = feature_dict_meta[sample][geneID][utr_length_to_include:-utr_length_to_include]
            sample_reads = sample_cds_reads[window[0]: window[1]]

            rna_cds_reads = feature_dict_meta[rn][geneID][utr_length_to_include:-utr_length_to_include]
            rn_reads = rna_cds_reads[window[0]: window[1]]

            geneLength = location_dict[geneID][1] - location_dict[geneID][0]
            
            # calculate rpkm at each position
            rpkm_sample = [i/geneLength/int(total_read_dict[sample])*10**9 for i in sample_reads]
            #rpkm_rna = [i/geneLength/int(total_read_dict[sample])*10**9 for i in rna_reads]

            #calculate TE at each position (rpkm sample/ rpkm rnaseq sample)
            temp_array.append([i/geneLength/int(total_read_dict[sample])*10**9 for i in sample_reads])
    
            # y axis values    
            # since temp_array is a single list, np.mean on axis 0 gives an array of the same values but with appropriate dimensions
            y_vals = np.mean(np.array(temp_array), axis=0)
            
            # calculate TE by normalizing to RPKM of gene (or RPKM of utr region)
            for j in range(0, len(y_vals)):
                y_vals[j] = y_vals[j] / (sum(rna_cds_reads) / geneLength/int(total_read_dict[sample])*10**9)
           

            plotting_dict[sample] = y_vals

    # list of x axis values, same length as reads
    x_vals = np.arange(window[0], window[1], 1)
    for sample in samples_to_consider:
        if sample in user_inputs['rnaSeq_pairs'].keys():
            ax.plot(x_vals, plotting_dict[sample], alpha=0.8, label=sample)

    
    seq = genome_dict['full_sequence_dict'][geneID][utr_length_to_include:-utr_length_to_include]
    aaSeq = str(Seq.translate(seq))

    # label with aa & tri-nt sequence at every 3rd nt
    x_ticks = [x for x in range(window[0], window[1]) if x % 3==0]
    x_labels = []
    for i in range(window[0], window[1]):
        if i % 3 == 0:
            x_labels.append(str(int(i/3)+1) +  '\n' + aaSeq[int(i/3)] + '\n' + str(seq[i:i+3])) 
    
    ax.legend(fontsize=16, loc='center left', bbox_to_anchor=(1, 0.5))
    
    ax.set_ylabel('TE (rpkm RiboSeq / rpkm RNAseq)', fontsize=16)    
    ax.set_xlabel(geneID.split('_')[2], style = 'italic', fontsize=16)
    ax.set_title('{} (aa {}-{})'.format(geneID.split('_')[2], int(window[0]/3)+1, int(window[1]/3)), fontsize=20)
    ax.set_xticks(x_ticks, x_labels)
    plt.savefig('{}/{}_TE_aa-{}-{}.pdf'.format(outdir, geneID.split('_')[2], int(window[0]/3)+1, int(window[1]/3)), bbox_inches='tight')

    

def whole_gene_window_rpm(geneID, samples_to_consider, total_read_dict, feature_dict_meta, genome_dict, outdir, user_inputs):

    print('Plot whole gene - RPM normalized - %s' % geneID)

    utr_length_to_include = user_inputs['utr_length_to_include']
    location_dict = genome_dict['location_dict']
    plotting_dict = {}

    fig,ax = plt.subplots()
    
    for sample in samples_to_consider:
        
        temp_array = []
    
        reads = feature_dict_meta[sample][geneID]
        
        # calculate RPM at each position
        temp_array.append([i/int(total_read_dict[sample])*10**6 for i in reads])
        
        # y axis values    
        # since temp_array is a single list, np.mean on axis 0 gives an array of the same values but with appropriate dimensions
        y_vals = np.mean(np.array(temp_array), axis=0)
        plotting_dict[sample] = y_vals
        
    
    # list of x axis values
    x_vals = np.arange(-utr_length_to_include, len(reads)-utr_length_to_include, 1)
       
    for sample in samples_to_consider:
        
        ax.plot(x_vals, plotting_dict[sample], alpha=0.8, label=sample)
    

    cds_seq = genome_dict['full_sequence_dict'][geneID][utr_length_to_include:-1*utr_length_to_include]
    seq = genome_dict['full_sequence_dict'][geneID]
    
    ax.legend(fontsize=16, loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_ylabel('reads per million (RPM)', fontsize=16)    
    ax.set_xlabel(geneID.split('_')[2], style = 'italic', fontsize=16)
    ax.set_title(geneID.split('_')[2], fontsize=20)
    
    # set x-axis ticks/labels at every 30 nt
    ax.xaxis.set_major_locator(MultipleLocator(30))
    lab = [int(i) for i in ax.get_xticks()]
    ax.set_xticks(ax.get_xticks(), lab, rotation=45)

    # bar & cutoffs to indicate CDS region
    min_y = min([int(i) for i in ax.get_yticks()])
    ax.annotate('', xy=(len(cds_seq), -0.15), xycoords=('data', 'axes fraction'),xytext=(0, -0.15), 
                arrowprops=dict(facecolor='grey', edgecolor='black',alpha=0.2, width=15, headwidth=20),annotation_clip=False)

    ax.axvline(len(cds_seq), c='r', linestyle ='--', linewidth=1)
    ax.axvline(0,  c='r', linestyle ='--', linewidth=1)
    
    plt.savefig('{}/{}_rpm_gene.pdf'.format(outdir, geneID.split('_')[2]), bbox_inches='tight')



def whole_gene_window_rpkm(geneID, samples_to_consider, total_read_dict, feature_dict_meta, genome_dict, outdir, user_inputs):

    print('Plot whole gene - RPKM normalized - %s' % geneID)

    utr_length_to_include = user_inputs['utr_length_to_include']
    location_dict = genome_dict['location_dict']
    plotting_dict = {}

    fig,ax = plt.subplots()
    
    for sample in samples_to_consider:
        
        temp_array = []
    
        reads = feature_dict_meta[sample][geneID]
        geneLength = location_dict[geneID][1] - location_dict[geneID][0]
        
        # calculate RPKM at each position
        temp_array.append([i/geneLength/int(total_read_dict[sample])*10**9 for i in reads])
        
        # y axis values    
        # since temp_array is a single list, np.mean on axis 0 gives an array of the same values but with appropriate dimensions
        y_vals = np.mean(np.array(temp_array), axis=0)
        plotting_dict[sample] = y_vals
        
    
    # list of x axis values
    x_vals = np.arange(-utr_length_to_include, len(reads)-utr_length_to_include, 1)
       
    for sample in samples_to_consider:
        
        ax.plot(x_vals, plotting_dict[sample], alpha=0.8, label=sample)
    

    cds_seq = genome_dict['full_sequence_dict'][geneID][utr_length_to_include:-1*utr_length_to_include]
    seq = genome_dict['full_sequence_dict'][geneID]
    
    ax.legend(fontsize=16, loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_ylabel('reads per kilobase per million (RPKM)', fontsize=16)    
    ax.set_xlabel(geneID.split('_')[2], style = 'italic', fontsize=16)
    ax.set_title(geneID.split('_')[2], fontsize=20)
    
    # set x-axis ticks/labels at every 30 nt
    ax.xaxis.set_major_locator(MultipleLocator(30))
    lab = [int(i) for i in ax.get_xticks()]
    ax.set_xticks(ax.get_xticks(), lab, rotation=45)

    # bar & cutoffs to indicate CDS region
    min_y = min([int(i) for i in ax.get_yticks()])
    ax.annotate('', xy=(len(cds_seq), -0.15), xycoords=('data', 'axes fraction'),xytext=(0, -0.15), 
                arrowprops=dict(facecolor='grey', edgecolor='black',alpha=0.2, width=15, headwidth=20),annotation_clip=False)

    ax.axvline(len(cds_seq), c='r', linestyle ='--', linewidth=1)
    ax.axvline(0,  c='r', linestyle ='--', linewidth=1)
    
    plt.savefig('{}/{}_rpkm_gene.pdf'.format(outdir, geneID.split('_')[2]), bbox_inches='tight')




def whole_gene_window_TE(geneID, samples_to_consider, total_read_dict, feature_dict_meta, genome_dict, outdir, user_inputs,RNAseq=False):

    if RNAseq == False:
        return(whole_gene_window_rpm(geneID, samples_to_consider, total_read_dict, feature_dict_meta, genome_dict, outdir, user_inputs))
    else:
        print('Plot whole gene - RPKM & TE normalized - %s' % geneID)

    utr_length_to_include = user_inputs['utr_length_to_include']
    location_dict = genome_dict['location_dict']
    plotting_dict = {}

    
    fig,ax = plt.subplots()
    
    for sample in samples_to_consider:
        if sample in user_inputs['rnaSeq_pairs'].keys():
            rn = user_inputs['rnaSeq_pairs'][sample]
            
            temp_array = []
        
            sample_reads = feature_dict_meta[sample][geneID]
            rn_reads = feature_dict_meta[rn][geneID]
            geneLength = location_dict[geneID][1] - location_dict[geneID][0]
            
            # calculate RPKM at each position
            rpkm_sample = [i/geneLength/int(total_read_dict[sample])*10**9 for i in sample_reads]

            # calculate TE at each position (rpkm sample/ rpkm rnaseq sample)
            temp_array.append(rpkm_sample)
            
            # y axis values    
            # since temp_array is a single list, np.mean on axis 0 gives an array of the same values but with appropriate dimensions
            y_vals = np.mean(np.array(temp_array), axis=0)

             # calculate TE by normalizing to RPKM of gene (or RPKM of utr region)
            for j in range(0, utr_length_to_include):
                y_vals[j] = y_vals[j] / (sum(rn_reads[0:utr_length_to_include]) / utr_length_to_include/int(total_read_dict[sample])*10**9)
            for j in range(utr_length_to_include, geneLength + utr_length_to_include):
                y_vals[j] = y_vals[j] / (sum(rn_reads[utr_length_to_include:geneLength + utr_length_to_include]) / geneLength/int(total_read_dict[sample])*10**9)
            for j in range(geneLength + utr_length_to_include, geneLength + 2*utr_length_to_include):
                y_vals[j] = y_vals[j] / (sum(rn_reads[geneLength + utr_length_to_include: geneLength + 2*utr_length_to_include]) / utr_length_to_include/int(total_read_dict[sample])*10**9)
            
            plotting_dict[sample] = y_vals
            
        
    # list of x axis values
    x_vals = np.arange(-utr_length_to_include, len(sample_reads)-utr_length_to_include, 1)
       
    for sample in samples_to_consider:
        if sample in user_inputs['rnaSeq_pairs'].keys():
        
            ax.plot(x_vals, plotting_dict[sample], alpha=0.8, label=sample)
    

    cds_seq = genome_dict['full_sequence_dict'][geneID][utr_length_to_include:-1*utr_length_to_include]
    seq = genome_dict['full_sequence_dict'][geneID]
    
    ax.legend(fontsize=16, loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_ylabel('TE (RPKM sample / RPKM RNAseq)', fontsize=16)    
    ax.set_xlabel(geneID.split('_')[2], style = 'italic', fontsize=16)
    ax.set_title(geneID.split('_')[2], fontsize=20)
    
    # set x-axis ticks/labels at every 30 nt
    ax.xaxis.set_major_locator(MultipleLocator(30))
    lab = [int(i) for i in ax.get_xticks()]
    ax.set_xticks(ax.get_xticks(), lab, rotation=45)

    # bar & cutoffs to indicate CDS region
    min_y = min([int(i) for i in ax.get_yticks()])
    ax.annotate('', xy=(len(cds_seq), -0.15), xycoords=('data', 'axes fraction'),xytext=(0, -0.15), 
                arrowprops=dict(facecolor='grey', edgecolor='black',alpha=0.2, width=15, headwidth=20),annotation_clip=False)

    ax.axvline(len(cds_seq), c='r', linestyle ='--', linewidth=1)
    ax.axvline(0,  c='r', linestyle ='--', linewidth=1)
    plt.savefig('{}/{}_TE_gene.pdf'.format(outdir, geneID.split('_')[2]), bbox_inches='tight')
    

def stack_plot(geneID, conditions_to_consider, outdir, feature_dict_meta, total_read_dict, genome_dict, user_inputs):
    
    geneID = analysis.gene_lookup(geneID, genome_dict)
    location_dict = genome_dict['location_dict']
    geneLength = location_dict[geneID][1] - location_dict[geneID][0]
    utr_length_to_include = user_inputs['utr_length_to_include']
    
    fig, axes = plt.subplots(nrows = len(conditions_to_consider), sharey=True, dpi=600)
    fig.supylabel('Average RPKM')
    fig.supxlabel('Gene Length (nt)')
    fig.suptitle(geneID.split('_')[-1])

    # list of x axis values
    x_vals = np.arange(-utr_length_to_include, geneLength + utr_length_to_include, 1)

    for i in range(0, len(conditions_to_consider)):
        ax = axes[i]
        condition = conditions_to_consider[i]

        RPKM_array = []

        for sample in user_inputs['condition_info'][condition]:
            temp_array = []
            
            reads = feature_dict_meta[sample][geneID]

            # calculate RPKM at each position
            temp_array.append([i/geneLength/int(total_read_dict[sample])*10**9 for i in reads])
  
            # since temp_array is a single list, np.mean on axis 0 gives an array of the same values but with appropriate dimensions
            sample_vals = np.mean(np.array(temp_array), axis=0)

            
            # y vals
            # calculate
            RPKM_array.append(sample_vals)
            avg_RPKM = np.mean(RPKM_array, axis = 0)
        
        ax.plot(x_vals, avg_RPKM, label = condition)
        ax.grid(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        if i != len(conditions_to_consider)-1:
            ax.set_xticks([])
            ax.axvline(0, ymin=-1.2, linestyle=':', linewidth = 1, c='black', clip_on=False)
            ax.axvline(geneLength, ymin=-1.2, linestyle=':', linewidth = 1, c='black', clip_on=False)
        else:
            ax.axvline(0, linestyle=':', linewidth = 1, c='black')
            ax.axvline(geneLength, linestyle=':', linewidth = 1, c='black')
        
    
        ax.annotate(condition, xy=(0.6,0.9),xycoords='axes fraction')

    plt.savefig('{}/{}_RPKM_stack.pdf'.format(outdir, geneID.split('_')[-1]), bbox_inches='tight')



def stack_plot_TE(geneID, conditions_to_consider, outdir, feature_dict_meta, total_read_dict, genome_dict, user_inputs, pos=()):

    geneID = analysis.gene_lookup(geneID, genome_dict)
    location_dict = genome_dict['location_dict']
    geneLength = location_dict[geneID][1] - location_dict[geneID][0]
    utr_length_to_include = user_inputs['utr_length_to_include']
    
    fig, axes = plt.subplots(nrows = len(conditions_to_consider), sharey=True, dpi=600)
    fig.supylabel('Average TE (RPKM Ribo / RPKM RNA)')
    fig.supxlabel('Gene Length (nt)')
    fig.suptitle(geneID.split('_')[-1])

    # list of x axis values
    x_vals = np.arange(-utr_length_to_include, geneLength + utr_length_to_include, 1)

    for i in range(0, len(conditions_to_consider)):
        ax = axes[i]
        condition = conditions_to_consider[i]

        TE_array = []

        for sample in user_inputs['condition_info'][condition]:
            
            temp_array = []
            
            reads = feature_dict_meta[sample][geneID]
            
            rn_sample = user_inputs['rnaSeq_pairs'][sample]
            rn_reads = feature_dict_meta[rn_sample][geneID]

            # calculate RPKM at each position
            temp_array.append([i/geneLength/int(total_read_dict[sample])*10**9 for i in reads])
  
            # since temp_array is a single list, np.mean on axis 0 gives an array of the same values but with appropriate dimensions
            sample_vals = np.mean(np.array(temp_array), axis=0)


            # calculate TE by normalizing to RPKM of gene (or RPKM of utr region)
            for j in range(0, utr_length_to_include):
                sample_vals[j] = sample_vals[j] / (sum(rn_reads[0:utr_length_to_include]) / utr_length_to_include/int(total_read_dict[sample])*10**9)
            for j in range(utr_length_to_include, geneLength + utr_length_to_include):
                sample_vals[j] = sample_vals[j] / (sum(rn_reads[utr_length_to_include:geneLength + utr_length_to_include]) / geneLength/int(total_read_dict[sample])*10**9)
            for j in range(geneLength + utr_length_to_include, geneLength + 2*utr_length_to_include):
                sample_vals[j] = sample_vals[j] / (sum(rn_reads[geneLength + utr_length_to_include: geneLength + 2*utr_length_to_include]) / utr_length_to_include/int(total_read_dict[sample])*10**9)


            # y vals
            # calculate
            TE_array.append(sample_vals)
            avg_TE = np.mean(TE_array, axis = 0)
        
        ax.plot(x_vals, avg_TE, label = condition)
        ax.grid(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
                   
        if i != len(conditions_to_consider)-1:
            ax.set_xticks([])
            #ax.set_ylim(0, 0.2)
            ax.axvline(0, ymin=-1.2, linestyle=':', linewidth = 1, c='black', clip_on=False)
            ax.axvline(geneLength, ymin=-1.2, linestyle=':', linewidth = 1, c='black', clip_on=False)
        else:
            ax.axvline(0, linestyle=':', linewidth = 1, c='black')
            ax.axvline(geneLength, linestyle=':', linewidth = 1, c='black')

        if len(pos)!=0:
            ax.axvline(pos[0], linestyle=':', linewidth = 1, c='red')
            ax.axvline(pos[1], linestyle=':', linewidth = 1, c='red')

    
        ax.annotate(condition, xy=(0.6,0.9),xycoords='axes fraction')
    

    plt.savefig('{}/{}_TE_stack.pdf'.format(outdir, geneID.split('_')[-1]), bbox_inches='tight')



###########################################################################################################################################
########################################## GENE PROFILING #################################################################################
###########################################################################################################################################







def gene_intraConditionCorrelation(df_master, user_inputs, outdir, RNAseq=False):
    #just using the first 2 samples per condition

    if RNAseq == False:
        datatype = 'RPKM'
    else:
        datatype = 'TE'
    
    
    for condition in user_inputs['condition_info']:
        if user_inputs['condition_type'][condition] == 'RiboSeq':
            compare_x = user_inputs['condition_info'][condition][0] + '_%s' % datatype
            compare_y = user_inputs['condition_info'][condition][1] + '_%s' % datatype
            temp_df = df_master[df_master[[compare_x , compare_y]].isnull().any(axis=1)==False]
        
            rho, p = stats.spearmanr(temp_df[compare_x], temp_df[compare_y])
            n = len(temp_df.index)
        
            fig, ax = plt.subplots()
            ax.loglog(temp_df[compare_x], temp_df[compare_y], c='k', alpha=0.2, marker='o', linestyle='')
            ax.set_xlabel(compare_x)
            ax.set_ylabel(compare_y)
            ax.set_title('Correlation: {:.3f}, (n= {})'.format(rho, n))
            plt.savefig('{}/{}_vs_{}_comparison.pdf'.format(outdir, compare_x, compare_y ), bbox_inches='tight')



def gene_interConditionCorrelation(df_master, user_inputs, outdir, RNAseq=False):
    
    if RNAseq == False:
        datatype = 'RPKM'
    else:
        datatype = 'TE'
        
    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    expConditions = [condition  for condition in user_inputs['condition_info'] if (condition != WT and user_inputs['condition_type'][condition] =='RiboSeq')]

    compare_WT = WT + '_avg%s' % datatype

    for condition in expConditions:
        compare_condition = condition + '_avg%s' % datatype
        
        temp_df = df_master[df_master[[compare_condition, compare_WT ]].isnull().any(axis=1)==False]
        rho, p = stats.spearmanr(temp_df[compare_condition], temp_df[compare_WT] )
        n = len(temp_df.index)
        
        fig, ax = plt.subplots()
        ax.loglog(temp_df[compare_condition], temp_df[compare_WT], c='k', alpha=0.2, marker='o', linestyle='')

        '''
        # in order to annotate specific genes, activate/modify this chunk
        
        for gene in ['typA', 'yjjK', 'dcuA', 'ispH', 'rhlB', 'hns']:

            xy = [(temp_df[compare_condition][k], temp_df[compare_WT][k]) for k in temp_df[compare_condition].keys() if gene in k][0]
            ax.loglog(xy[0], xy[1],'ro',alpha=0.5, label=gene)
            #ax.text(xy[0] * 1.1, y[i] * 1.1, label, fontsize=12)
            ax.annotate(
                    gene,                # Text to display
                    xy,      # The point to annotate (x, y)
                    xytext=(xy[0] * 0.5, xy[1] * 10), # Position of the text
                    arrowprops=dict(            # Arrow properties
                        color='red', 
                        arrowstyle='->', 
                        connectionstyle='arc3'
                    )
                )
        '''
            
        ax.set_ylabel(compare_WT)
        ax.set_xlabel(compare_condition)
        ax.set_title('Correlation: {:.3f}, (n= {})'.format(rho, n))
        plt.savefig('{}/{}_vs_{}_comparison.pdf'.format(outdir, compare_condition, compare_WT), bbox_inches='tight')

def gene_histFC(df_master, user_inputs, outdir, RNAseq=False):
    
    if RNAseq == False:
        datatype = 'RPKM'
    else:
        datatype = 'TE'

    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    expConditions = [condition  for condition in user_inputs['condition_info'] if (condition != WT and user_inputs['condition_type'][condition] =='RiboSeq')]
    
    for condition in expConditions:
        
        fig, ax = plt.subplots()
        plt.hist(df_master['{}_{}_log2FC'.format(condition, datatype)], bins = 30, facecolor='gainsboro', edgecolor='dimgray')
        plt.xlabel('Fold Change %s vs %s - '  % (condition, WT) \
                   + '$log2(\ \dfrac{ %s_{%s} }{ %s_{%s} }\ )$' % (datatype, condition, datatype, WT))
        plt.ylabel('Number of Genes')
        plt.axvline(np.nanpercentile(df_master['{}_{}_log2FC'.format(condition, datatype)], 5) , c='r', linestyle = ':')
        plt.axvline(np.nanpercentile(df_master['{}_{}_log2FC'.format(condition, datatype)], 95) , c='r', linestyle = ':')
        plt.savefig('{}/{}-{}_{}_log2FC_hist.pdf'.format(outdir, condition, WT, datatype), bbox_inches='tight')






###########################################################################################################################################
########################################## CODON PROFILING ################################################################################
###########################################################################################################################################




def codon_waterfall(sorted_log2FC_pauseScore, sensitive_dict, resistant_dict, outdir, user_inputs):

    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    
    for condition in sorted_log2FC_pauseScore:
        fig, ax = plt.subplots(figsize=(6, 7))
    
        lists = sorted_log2FC_pauseScore[condition].items()
        x, y = zip(*lists) # unpack a list of pairs into two tuples
    
        plt.scatter(range(0, len(x)), y, color = 'k', s=2.5, linewidths=0)  ## dot for every individual codon
        #plt.plot(range(0, len(x)), y, color = 'k')   ## makes a smooth line
    
        plt.title(condition)
        plt.ylabel('log2 ( PauseScore %s / PauseScore %s )' % (condition, WT))
        plt.text(len(x)*.95, 0.1,'codons', horizontalalignment='left')
        
        
        ax.xaxis.set_major_locator(MultipleLocator(40000))
        ax.spines[['right', 'top']].set_visible(False)
        ax.spines['bottom'].set_position('zero')

        xmin, xmax =ax.get_xlim()
        bar_ratio = len(x) / len(range(int(xmin), int(xmax)))

        norm = colors.Normalize(max(sorted_log2FC_pauseScore[condition].values()), min(sorted_log2FC_pauseScore[condition].values()))
        cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap = BuViPi), ax=ax, orientation='horizontal', pad=0.01, shrink=bar_ratio, label = 'Stall Strength')
        cbar.set_ticks([])
       
        ax.grid(False)
        plt.savefig(outdir+'%s-%s_log2FC_codonWaterfall.pdf' % (condition, WT), bbox_inches='tight')

        # add lines for 5 and 95 percentile
        ax.axvline(len(x)-len(resistant_dict[condition]), color ='r', linestyle = ':')
        ax.axvline(len(sensitive_dict[condition]), color ='r', linestyle = ':')

        plt.savefig(outdir+'%s-%s_log2FC_codonWaterfall+cutoffs.pdf' % (condition, WT), bbox_inches='tight')


def codon_histFC(sorted_log2FC_pauseScore, outdir, user_inputs):
    
    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    
    for condition in sorted_log2FC_pauseScore:
        fig, ax = plt.subplots()
        plt.hist(sorted_log2FC_pauseScore[condition].values(), bins=20, facecolor='gainsboro', edgecolor='dimgray')
        plt.xlabel('log2 ( PauseScore %s / PauseScore %s )' % (condition, WT))
        plt.ylabel('Number of Codons')
        plt.title(condition)
        plt.savefig(outdir+'%s-%s_log2FC_codonHist.pdf' % (condition, WT), bbox_inches='tight')

        # add lines for 5 and 95 percentile
        ax.axvline(np.percentile(list(sorted_log2FC_pauseScore[condition].values()), 95), color ='r', linestyle = ':')
        ax.axvline(np.percentile(list(sorted_log2FC_pauseScore[condition].values()), 5), color ='r', linestyle = ':')

        plt.savefig(outdir+'%s-%s_log2FC_codonHist+cutoffs.pdf' % (condition, WT), bbox_inches='tight')

        

colors2 = ['#fe6100', '#dc267e', '#795ef0', '#56b4e9', '#009e74']

def codon_frequency(codon_freq, outdir, user_inputs):
    
    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    
    i= 0
    codon_color_map=[]
    for aa in sorted(aaTable):
        for codon in aaTable[aa]:
            codon_color_map.append( colors2[i%len(colors2)] )
        i+=1
    
    
    for condition in codon_freq['P']:
        maxy= max(list(codon_freq['E'][condition].values()) + list( codon_freq['P'][condition].values()) + list( codon_freq['A'][condition].values()))
        miny= min(list(codon_freq['E'][condition].values()) + list( codon_freq['P'][condition].values()) + list( codon_freq['A'][condition].values()))
        
        for position in codon_freq:
            
            fig, ax = plt.subplots( figsize=(20, 5))
    
            sorted_lists = sorted(codon_freq[position][condition].items(), key=lambda x: gencode[x[0]])[3:]
            codons, vals = zip(*sorted_lists) # unpack a list of pairs into two tuples
            
            x = np.arange(len(codons))
    
            plt.bar(x, vals, color=codon_color_map[3:], label=codons)
            plt.title(condition + '  -  ' + position + ' Site Codon Frequency',  fontsize=16)
            plt.ylabel('avg log2FC vs %s'% WT)
            ax.set_xticks(x, codons, rotation = 90)
            ycount = len(ax.get_yticks())
            plt.ylim(-1.5, 1.5)
            ax.set_yticks(np.arange(-1.5, 1.6, 0.25),np.arange(-1.5, 1.6, 0.25)) 
            
            
            
            ax.grid(False, axis='x')
            #ax.tick_params(axis='x', colors=codon_color_map[3:])
           
            # lines between the AAs:
            sec2 = ax.secondary_xaxis(location=0)
            divisions = [-0.5]
            for aa in aaTable :
                if aa != '*':
                    divisions.append(divisions[-1] + len(aaTable[aa]))
            sec2.set_xticks(divisions, labels=[])
            sec2.tick_params('x', length=50, width=1.5)
    
            # label the AAs:
            sec = ax.secondary_xaxis(location=0)
            sec.set_xticks([np.mean(divisions[i:i+2]) for i in range(len(divisions)-1)], labels= aaList[1:])
            sec.tick_params('x', length=40, width=0)
    
            plt.savefig(outdir+'%s-%s_log2FC_%s-site_codonFreq.pdf' % (condition, WT, position), bbox_inches='tight')    



def aa_frequency(aa_freq, outdir, user_inputs):
    
    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    
    aa_color_map=[]
    for i, aa in enumerate(sorted(aaTable)):
        aa_color_map.append( colors2[i%len(colors2)] )
        i+=1
    
    
    for condition in aa_freq['P']:
        maxy= max(list(aa_freq['E'][condition].values()) + list( aa_freq['P'][condition].values()) + list( aa_freq['A'][condition].values()))
        miny= min(list(aa_freq['E'][condition].values()) + list( aa_freq['P'][condition].values()) + list( aa_freq['A'][condition].values()))
        
        for site in aa_freq:
            
            fig, ax = plt.subplots()
    
            sorted_lists = sorted(aa_freq[site][condition].items())[1:]
            aas, vals = zip(*sorted_lists) # unpack a list of pairs into two tuples
            
            x = np.arange(len(aas))
    
            plt.bar(x, vals, color=aa_color_map, label=aas)
            plt.title(condition + '  -  ' + site + ' Site Amino Acid Frequency',  fontsize=16)
            plt.ylabel('avg log2FC vs %s'% WT)
            ax.set_xticks(x, aas)
            ycount = len(ax.get_yticks())
            plt.ylim(-1.5, 1.5)
            ax.set_yticks(np.arange(-1.5, 1.6, 0.25),np.arange(-1.5, 1.6, 0.25)) 
            ax.grid(False, axis='x')
    
            plt.savefig(outdir+'%s-%s_log2FC_%s-site_aaFreq.pdf' % (condition, WT, site), bbox_inches='tight')    



#######################################################################################################################################################

def seq_hmap(codon_cpm_meta, genome_dict, outdir, user_inputs): 

    positions = range(-1,-12, -1)
    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    aaList = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q',
              'C', 'G', 'P', 'A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V', '*']

    print('heatmap by positions: ', [x+2 for x in positions][::-1])
    print('excluding stop codons except at A site')
    
    seqPos_dict = {}
    for condition in codon_cpm_meta:
        seqPos_dict[condition] = {}
        for sample in codon_cpm_meta[condition]:
            seqPos_dict[condition][sample] = {x:{y:[] for y in aaList} for x in positions}
            for codon_id in codon_cpm_meta[condition][sample]:
                seq = analysis.get_peak_seq_single(codon_id, genome_dict, user_inputs)
                cpm = codon_cpm_meta[condition][sample][codon_id]

                # exclude stop codon seq?
                if '*' not in seq[:-1]:
        
                    #make a list of cpm for each aa at positions -11 - +1 (for each sample individually)
                    for i, pos in enumerate(positions):
                        if i< len(seq):
                            #print(seq, pos, seq[pos])
                            seqPos_dict[condition][sample][pos][seq[pos]].append(cpm)
        
            # calc average cpm for each aa (for each sample individually)
            for pos in seqPos_dict[condition][sample]:
                for aa in seqPos_dict[condition][sample][pos]:
                    aaAvg = np.mean(seqPos_dict[condition][sample][pos][aa])
                    seqPos_dict[condition][sample][pos][aa] = aaAvg

    
    for condition in seqPos_dict: 
        for C in [condition, WT]:
            seqPos_dict[condition][C + '_avg'] = {x:{y:[] for y in aaList} for x in positions}
            seqPos_dict[condition][C + '_log2'] = {x:{y:0 for y in aaList} for x in positions}
        seqPos_dict[condition][condition + '_log2FC'] = {x:{y:0 for y in aaList} for x in positions}
        for pos in positions:
            for aa in aaList:
                for sample in user_inputs['condition_info'][condition]:
                    seqPos_dict[condition][condition + '_avg'][pos][aa].append(seqPos_dict[condition][sample][pos][aa])
                for sample in user_inputs['condition_info'][WT]:
                    seqPos_dict[condition][WT + '_avg'][pos][aa].append(seqPos_dict[condition][sample][pos][aa])
                    
                for C in [condition, WT]:
                    
                    # calc average cpm for each aa/position (avg between replicates)
                    avg = np.mean(seqPos_dict[condition][C + '_avg'][pos][aa])
                    seqPos_dict[condition][C + '_avg'][pos][aa] = avg
    
                    #calc log2 fold change compared to DMSO (from condition averages)
                    log2 = np.log2(avg)
                    seqPos_dict[condition][C + '_log2'][pos][aa] = log2
    
                seqPos_dict[condition][condition + '_log2FC'][pos][aa] = seqPos_dict[condition][condition + '_log2'][pos][aa] - seqPos_dict[condition][WT + '_log2'][pos][aa]
                
    for condition in seqPos_dict:
        df = pd.DataFrame(seqPos_dict[condition][condition + '_log2FC'])
        df = df.rename(columns = {x: x+2 for x in positions})
        df = df[df.columns[::-1]]
        print('\n', condition) 
        display(df)
        
        ax = sns.heatmap(df, cmap='coolwarm') #, vmin=-1.25, vmax=1.25)
        ax.set_title('Log2FC RPM %s vs RPM %s' %(condition, WT))
        ax.grid(False)
        for _, spine in ax.spines.items():
            spine.set_visible(True)

        plt.savefig(outdir+'%s-%s_log2FC_seqHmap.pdf' % (condition, WT), bbox_inches='tight') 
        plt.show()
        plt.close()



def pairwise_hmap2(condition, pos, codon_cpm_meta, outdir, genome_dict, user_inputs, ax=None, grid=False, transform='np.log2'):
    sites = {'A': 1, 'P': 0, 'E': -1}
    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    
    # "if x position, then y position"
    pairwise_dict = {}
    xpos = [sites[pos[0]]-2 if pos[0] in sites else int(pos[0])-2][0]
    ypos = [sites[pos[1]]-2 if pos[1] in sites else int(pos[1])-2][0]
    #print(pos,xpos, ypos)

    for sample in codon_cpm_meta[condition]:
        pairwise_dict[sample] = {x: {y:[] for y in aaList} for x in aaList}
        for codon_id in codon_cpm_meta[condition][sample]:
            seq = analysis.get_peak_seq_single(codon_id, genome_dict, user_inputs)
            cpm = codon_cpm_meta[condition][sample][codon_id]

            # exclude stop codon seq?
            if '*' not in seq[:-1]:
                if ypos >= -len(seq) and xpos >= -len(seq):
                    # list cpm by position x AND position y (for each sample individually)
                    pairwise_dict[sample][seq[xpos]][seq[ypos]].append(cpm)

        #average cpm for each aa pair (for each sample individually)
        for aaX in pairwise_dict[sample]:
            for aaY in pairwise_dict[sample][aaX]:
                pairAvg = np.mean(pairwise_dict[sample][aaX][aaY])
                pairwise_dict[sample][aaX][aaY] = pairAvg
    
    for C in [condition, WT]:
        pairwise_dict[C + '_avg']= {x:{y:[] for y in aaList} for x in aaList}
        #print(condition, pairwise_dict.keys())
        for x in aaList:
            for y in aaList:
                for sample in user_inputs['condition_info'][C]:
                    pairwise_dict[C + '_avg'][x][y].append(pairwise_dict[sample][x][y])

                # calc avg cpm for each aa pair (avg between replicates)
                avg = np.mean(pairwise_dict[C + '_avg'][x][y])
                pairwise_dict[C + '_avg'][x][y] = avg

    # calc fold change or log2 fold change compared to DMSO (from condition averages)
    if transform == None:
        pairwise_dict[condition + '_FC'] = {x:{y:[] for y in aaList} for x in aaList} 
        for x in aaList:
            for y in aaList:
                pairwise_dict[condition + '_FC'][x][y] = pairwise_dict[condition + '_avg'][x][y] / pairwise_dict[WT + '_avg'][x][y]
        
        df = pd.DataFrame(pairwise_dict[condition + '_FC'])
        #print('\n', condition) 
        #display(df)
        
        ax = sns.heatmap(df, ax=ax, cmap='coolwarm', xticklabels=True, yticklabels=True, square=True, cbar_kws= {'shrink': 0.5})
        
        ax.grid(False)
    
        for a in ax.figure.axes:
            if a.get_label() == '<colorbar>':
                a.set_ylabel('FC RPM %s vs RPM %s' % (condition, WT))
        for s in ax.spines.values():
            s.set_visible(True)

        if grid==False:
            ax.set_title('FC RPM %s vs RPM %s' %(condition, WT))
            plt.xlabel(pos[0] + ' site')
            plt.ylabel(pos[1] + ' site')
            plt.savefig(outdir+'%s-%s_%s:%s_FChmap.pdf' % (condition, WT, pos[1], pos[0]), bbox_inches='tight')
            plt.show()
            plt.close()
            
    elif transform == 'np.log2':
        pairwise_dict[condition + '_log2FC'] = {x:{y:[] for y in aaList} for x in aaList} 
        for x in aaList:
            for y in aaList:
                pairwise_dict[condition + '_log2FC'][x][y] = ( 
                    np.log2(pairwise_dict[condition + '_avg'][x][y]) 
                    - np.log2(pairwise_dict[WT + '_avg'][x][y]) )
                #pairwise_dict[condition + '_avg'][x][y] / pairwise_dict[WT + '_avg'][x][y]
        
        df = pd.DataFrame(pairwise_dict[condition + '_log2FC'])
        #print('\n', condition) 
        #display(df)
        
        ax = sns.heatmap(df, ax=ax, cmap='coolwarm', xticklabels=True, yticklabels=True, square=True, cbar_kws= {'shrink': 0.5})

        ax.grid(False)
    
        for a in ax.figure.axes:
            if a.get_label() == '<colorbar>':
                a.set_ylabel('log2FC RPM %s vs RPM %s' % (condition, WT))
        for s in ax.spines.values():
            s.set_visible(True)

        if grid==False:
            ax.set_title('Log2FC RPM %s vs RPM %s' %(condition, WT))
            plt.xlabel(pos[0] + ' site')
            plt.ylabel(pos[1] + ' site')
            plt.savefig(outdir+'%s-%s_%s:%s_log2FChmap.pdf' % (condition, WT, pos[1], pos[0]), bbox_inches='tight')
            plt.show()
            plt.close()

    return ax



def hmap_grid(condition, codon_cpm_meta, outdir, genome_dict, user_inputs, transform='np.log2'):
    pos = ['-2', 'E', 'P', 'A']
    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    grid_size = len(pos) - 1
    f, axes = plt.subplots(
        nrows=grid_size,
        ncols=grid_size,
        figsize=(grid_size * 5, grid_size * 5),
    )
    for i, c in enumerate(pos[:-1]):
        for j, r in enumerate(pos[1:]):
            ax = axes[i, j]
            if i > j:
                ax.set_visible(False)
                continue
            # if i == j:
            #    ax.set_ylabel(c)
            if i == 0:
                ax.set_title(r)
            pairwise_hmap2(condition, (r, c), codon_cpm_meta, outdir, genome_dict, user_inputs, ax, grid=True, transform=transform)
            ax.set_xlabel(r)
            ax.set_ylabel(c)
    print(condition)
    f.savefig(outdir+'%s-%s_%shmapGrid.pdf' % (condition, WT, ['' if transform==None else transform[3:]][0]), bbox_inches='tight')
    

def volcano_trimer(codon_counts_meta, outdir, genome_dict, user_inputs):
    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    trimerlist = []
    for aa1 in aaList:
        for aa2 in aaList:
            for aa3 in aaList:
                if aa1 != '*' and aa2 != '*':
                    trimerlist.append(aa1+aa2+aa3)
                    
    trimer_counts_dict = {}
    for condition in codon_counts_meta:
        trimer_counts_dict[condition] = {}
        for sample in codon_counts_meta[condition]:
            trimer_counts_dict[condition][sample] = {x: [] for x in trimerlist}
            for codon_id in codon_counts_meta[condition][sample]:
                seq = analysis.get_peak_seq_single(codon_id, genome_dict, user_inputs)
                counts = codon_counts_meta[condition][sample][codon_id]
    
               # exclude stop codon seq?
                if '*' not in seq[:-1]:
    
                    #make a list of cpm for each trimer at positions E-A-P (for each sample individually)
                    trimer_counts_dict[condition][sample][seq[-3:]].append(counts)
    
            # calc average cpm for each trimer (for each sample individually)
            for trimer in trimer_counts_dict[condition][sample]:
                trimerAvg = np.mean(trimer_counts_dict[condition][sample][trimer])
                trimer_counts_dict[condition][sample][trimer] = trimerAvg
    
    for condition in trimer_counts_dict:
        df = pd.DataFrame(trimer_counts_dict[condition])
        metadata ={'sample_id': [],'condition': [], 'group': [] }
        for sample in codon_counts_meta[condition].keys():
            metadata['sample_id'].append(sample)
            if sample in user_inputs['condition_info'][condition]:
                metadata['condition'].append(condition)
            elif sample in user_inputs['condition_info'][WT]:
                metadata['condition'].append(WT)
            for repnum in user_inputs['repNums']:
                if sample in user_inputs['repNums'][repnum]:
                    metadata['group'].append(repnum)
        metadata = pd.DataFrame(metadata)
        metadata = metadata.set_index('sample_id')
                
        dds = DeseqDataSet(counts=df.fillna(0)[codon_counts_meta[condition].keys()].T.astype(int), 
                           metadata=metadata, refit_cooks=True)
    
        dds.deseq2()
    
        stat_res = DeseqStats(dds, contrast=['condition', condition, WT])
        stat_res.summary()
    
        res = stat_res.results_df.sort_values(by='log2FoldChange', ascending=False)
        res['log10pvalue'] = -np.log10(res['pvalue'])
        res['log10padj'] = -np.log10(res['padj'])
        
        data = {'log2FC': res['log2FoldChange'], 'log10padj':res['log10padj']}
        df = pd.DataFrame(data)
        display(df)
        plt.figure(figsize=(10,10))
        plt.scatter(x=res['log2FoldChange'], y=res['log10padj'], s=10)
        plt.title('%s vs. %s: E:A' % (condition, WT))
        plt.xlabel('log2FoldChange')
        plt.ylabel('log10padj')
    
        
        for pos in df.index:
            if res['log10padj'][pos]>5 and res['log2FoldChange'][pos] > 2.5:
                plt.annotate(pos, xy =(res['log2FoldChange'][pos] + .025, res['log10padj'][pos] +.075))
            if res['log10padj'][pos]>5 and res['log2FoldChange'][pos] < -2.5:
                plt.annotate(pos, xy =(res['log2FoldChange'][pos] + .025, res['log10padj'][pos] +.075))
    
        #plt.show()
        plt.savefig(outdir+'%s-%s_volcano3.pdf' % (condition, WT), bbox_inches='tight')


def volcano_tetramer(codon_counts_meta, outdir, genome_dict, user_inputs):
    WT = ['DMSO' if 'DMSO' in user_inputs['condition_info'] else 'WT'][0]
    tetramerlist = []
    for aa1 in aaList:
        for aa2 in aaList:
            for aa3 in aaList:
                for aa4 in aaList:
                    if aa1 != '*' and aa2 != '*' and aa3!= '*':
                        tetramerlist.append(aa1+aa2+aa3+aa4)
    tetramer_counts_dict = {}
    for condition in codon_counts_meta:
        tetramer_counts_dict[condition] = {}
        for sample in codon_counts_meta[condition]:
            tetramer_counts_dict[condition][sample] = {x: [] for x in tetramerlist}
            for codon_id in codon_counts_meta[condition][sample]:
                seq = analysis.get_peak_seq_single(codon_id, genome_dict, user_inputs)
                counts = codon_counts_meta[condition][sample][codon_id]
    
               # exclude stop codon seq?
                if '*' not in seq[:-1]:
    
                    #make a list of cpm for each tetramer at positions -2-E-A-P (for each sample individually)
                    tetramer_counts_dict[condition][sample][seq[-4:]].append(counts)
    
            # calc average cpm for each tetramer (for each sample individually)
            for tetramer in tetramer_counts_dict[condition][sample]:
                tetramerAvg = np.mean(tetramer_counts_dict[condition][sample][tetramer])
                tetramer_counts_dict[condition][sample][tetramer] = tetramerAvg
    
    for condition in tetramer_counts_dict:
        df = pd.DataFrame(tetramer_counts_dict[condition])
        metadata ={'sample_id': [],'condition': [], 'group': [] }
        for sample in codon_counts_meta[condition].keys():
            metadata['sample_id'].append(sample)
            if sample in user_inputs['condition_info'][condition]:
                metadata['condition'].append(condition)
            elif sample in user_inputs['condition_info'][WT]:
                metadata['condition'].append(WT)
            for repnum in user_inputs['repNums']:
                if sample in user_inputs['repNums'][repnum]:
                    metadata['group'].append(repnum)
        metadata = pd.DataFrame(metadata)
        metadata = metadata.set_index('sample_id')
                
        #rounding averages to ints in order to use Deseq2
        dds = DeseqDataSet(counts=df.fillna(0)[codon_counts_meta[condition].keys()].T.astype(int), 
                           metadata=metadata,
                           refit_cooks=True)
    
        dds.deseq2()
    
        stat_res = DeseqStats(dds, contrast=['condition', condition, WT])
        stat_res.summary()
    
        res = stat_res.results_df.sort_values(by='log2FoldChange', ascending=False)
        res['log10pvalue'] = -np.log10(res['pvalue'])
        res['log10padj'] = -np.log10(res['padj'])
        
        data = {'log2FC': res['log2FoldChange'], 'log10padj':res['log10padj']}
        df = pd.DataFrame(data)
        display(df)
        plt.figure(figsize=(10,10))
        plt.scatter(x=res['log2FoldChange'], y=res['log10padj'], s=10)
        plt.title('%s vs. %s: -2:A' % (condition, WT))
        plt.xlabel('log2FoldChange')
        plt.ylabel('log10padj')
    
        
        for pos in df.index:
            if res['log10padj'][pos]>20 and res['log2FoldChange'][pos] > 2.5:
                plt.annotate(pos, xy =(res['log2FoldChange'][pos] + .025, res['log10padj'][pos] +.075))
    
        #plt.show()
        plt.savefig(outdir+'%s-%s_volcano4.pdf' % (condition, WT), bbox_inches='tight')
        