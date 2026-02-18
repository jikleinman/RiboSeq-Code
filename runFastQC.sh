#runFastQC.sh

#FastQC documentation @ https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/ 

#full path for text file containing relevant adapter names & sequences separated by tab
adapterfile=$1
#folder to write FastQC results
outdir=$2
#full path for fastq or fastq.gz file to be analyzed
infile=$3

fastqc -a $adapterfile -o $outdir $infile