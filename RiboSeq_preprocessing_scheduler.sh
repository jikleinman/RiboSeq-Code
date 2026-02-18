#!/bin/bash
#$ -S /bin/bash

#$ -o ~/JKDF08/wyntonlogs/JKDF08_preprocessing_log.txt
#$ -e ~/JKDF08/wyntonlogs/JKDF08_preprocessing_error.txt
#$ -cwd
#$ -pe smp 40
#$ -r y
#$ -j y
#$ -l mem_free=2G
#$ -l scratch=2G
#$ -l h_rt=12:00:00
#$ -N JKDF08_preprocessing


#export PATH=/salilab/diva1/home/anaconda/anaconda3/bin:$PATH
#source activate riboseq-env2
export PATH=/wynton/home/cbi/shared/software/CBI/miniconda3-23.5.2-0-py311/bin:$PATH
source activate riboseq-env

# the actual command here is:
#python RiboSeq_preprocessing.py filename1 filename2 filename3 ...
# but since i was having trouble with the environment following a wynton update, may need to include it in the python path - ex:
~/.conda/envs/riboseq-env/bin/python RiboSeq_preprocessing.py JKDF06_Rn2_DMSO_S75_L003_R1_001 JKDF06_Rn2_100xDZD_S77_L003_R1_001 JKDF06_Rn2_100xCZD_S79_L003_R1_001 JKDF06_Rn2_10xDZD_S76_L003_R1_001 JKDF06_Rn2_10xCZD_S78_L003_R1_001 JKDF06_Rn1_DMSO_S70_L003_R1_001 JKDF06_Rn1_100xDZD_S72_L003_R1_001 JKDF06_Rn1_100xCZD_S74_L003_R1_001 JKDF06_Rn1_10xDZD_S71_L003_R1_001 JKDF06_Rn1_10xCZD_S73_L003_R1_001 
 

#Rb1_DMSO Rb1_TZD-10 Rb1_TZD-50 Rb1_TZD-100 Rb2_DMSO Rb2_TZD-10 Rb2_TZD-50 Rb2_TZD-100


[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"