#RiboSeq_runCutadapt.sh

#cutadapt documentation @ https://cutadapt.readthedocs.io/en/stable/guide.html#basic-usage

#full path for input data fastq or fastq.gz
infile=$1
#full path for output (1) fastq.gz containing all passing reads
cutfile=$2
#full path for output (2) fastq.gz containing all reads below minLen
shortfile=$3
#full path for output (3) fastq.gz containing all reads above maxLen
longfile=$4
#full path for output (4) txt file containing cutadapt report
cutlog=$5
#3' adapter nucleotide sequence to be trimmed
threePrime=$6
#minimum read length AFTER trimming
minLen=$7
#maximum read length after trimming (also for untrimmed reads)
maxLen=$8
#if N>1: total allowable errors in adapter matching; if 0≤N<1: acceptable ratio of errors to adapter length 
acceptableError=$9
#number of cores to use in parallel processing - for .gz files, requires pigz to be installed (seems to be default on wynton)
numCores=${10}

#Trims the provided 3' adapter from all reads (both full & partial occurences, plus all downstream nucleotides)
#Where resulting trimmed read is not in the specified length range, saved to appropriate long/short file

cutadapt -a $threePrime -m $minLen -M $maxLen -j $numCores -e $acceptableError --too-short-output $shortfile --too-long-output $longfile -o $cutfile $infile >> $cutlog
