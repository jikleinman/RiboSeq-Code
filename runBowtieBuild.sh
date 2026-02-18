# runBowtieBuild.sh

# builds a Bowtie index from a set of DNA sequences; 
# outputs a set of 6 files with suffixes .1.ebwt, .2.ebwt, .3.ebwt, .4.ebwt, .rev.1.ebwt, and .rev.2.ebwt

# bowtie documentation: https://bowtie-bio.sourceforge.net/manual.shtml#the-bowtie-aligner

#reference genome, in fasta format
ref_in=$1

#name for output of the reference genome index
ebwt_base=$2

bowtie-build $ref_in $ebwt_base

