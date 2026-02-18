# RiboSeq_runBowtieAlign.sh

# v3 = 20240911 filters reads via a second fasta file to remove tRNA/rRNA/etc before whole genome alignment
# pair with Riboseq_preprocessing.ipynb/py w/ v5


#bowtie takes an index and a set of reads as input and outputs a list of alignments
#works best when aligning short reads to large genomes (e.g. human or mouse), though it supports arbitrarily small ref sequences and reads as long as 1024 bases
#documentation: https://bowtie-bio.sourceforge.net/manual.shtml#the-bowtie-aligner
#ebwt_base, infile, sam_outdir, bwt_outdir, alignSam, alignBwt, unalignFastq, alignFastq, maxFastq, bowtieLog, acceptableError, numCores

###################################################################################################################################################################

fname=$1
#path to initial fastq file being aligned
infile=$2
#outdir for initial (filtering) alignment
filtering_outdir=$3
#outdir for final (genome) alignment
genomeAlign_outdir=$4
#path to reference tRNA/rRNA index (-x)
filtering_ebwt=$5
#path to reference genome index (-x)
genome_ebwt=$6
#path to .txt file with timing and processing statistics
bowtieLog=$7
# -v option: Report alignments with at most <int> mismatches
# -v 2 enforces simpler end-to-end k-difference policy 
acceptableError=$8
# -p option: number of cores to use in parallel processing - must match scheduler
# --reorder option: ensures read order in SAM file matches order of appearance in infile when using parallel processing
numCores=$9

# -y option: Try as hard as possible to find valid alignments when they exist; slower than default

# -m option: Suppress all alignments for a particular read or pair if more than <int> reportable alignments exist for it
# for riboSeq -m = 1 in order to eliminate reads from rRNA 

# -t option: print timing statistics

# -a option: Report all valid alignments per read. If more than one valid alignment exists and the --best and --strata options are specified, then only those alignments belonging to the best alignment “stratum” will be reported.

# --best option: Make Bowtie guarantee that reported singleton alignments are “best” in terms of stratum 
#(i.e. number of mismatches, or mismatches in the seed in the case of -n mode) and in terms of the quality values at the mismatched position(s). 
#Stratum always trumps quality; e.g. a 1-mismatch alignment where the mismatched position has Phred quality 40 is preferred over a 2-mismatch alignment where the mismatched positions both have Phred quality 10

# --strata: If many valid alignments exist and are reportable (e.g. are not disallowed via the -k option) and they fall into more than one alignment “stratum”, report only those alignments that fall into the best stratum. By default, Bowtie reports all reportable alignments regardless of whether they fall into multiple strata. When --strata is specified, --best must also be specified.

# -S option: output will be in SAM format

echo 'input: '$infile' - filtering: '$(basename $filtering_ebwt)
echo $'\n'$fname' - filtering: '$(basename $filtering_ebwt) 2>&1 | tee -a $bowtieLog
(bowtie -v $acceptableError -t -y -m 1 -a --best --strata -x $filtering_ebwt $infile -S $filtering_outdir$fname'_filtering_match.sam' --un $filtering_outdir$fname'_filtering_removed.fastq' --al $filtering_outdir$fname'_filtering_match.fastq' --max $filtering_outdir$fname'_filtering_multi.fastq')  2>>$bowtieLog
#convert fastq outputs to .gz for storage
gzip $filtering_outdir$fname'_filtering_removed.fastq'
gzip $filtering_outdir$fname'_filtering_match.fastq'
gzip $filtering_outdir$fname'_filtering_multi.fastq'
echo filtered $fname

echo 'input: '$filtering_outdir$fname'_filtering_removed.fastq.gz'' - aligning to genome: '$(basename $genome_ebwt)
echo $'\n'$fname'_filtered - aligning to genome: '$(basename $genome_ebwt) 2>&1 | tee -a $bowtieLog
(bowtie -v $acceptableError -t -y -m 1 -a --best --strata -x $genome_ebwt $filtering_outdir$fname'_filtering_removed.fastq.gz' -S $genomeAlign_outdir$fname'_filtered_genome_match.sam' --un $genomeAlign_outdir$fname'_filtered_genome_nomatch.fastq' --al $genomeAlign_outdir$fname'_filtered_genome_match.fastq' --max $genomeAlign_outdir$fname'_filtered_genome_multi.fastq')  2>>$bowtieLog
#convert fastq outputs to .gz for storage
gzip $genomeAlign_outdir$fname'_filtered_genome_nomatch.fastq'
gzip $genomeAlign_outdir$fname'_filtered_genome_match.fastq'
gzip $genomeAlign_outdir$fname'_filtered_genome_multi.fastq'
echo aligned $fname


