#!/bin/bash

# RUN THIS SCRIPT IN ONE FOLDER OUTSIDE 01.RAW DATA (colias2023)

#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=TRIM
#################
#a file for job output, you can check job progress for each sample
#SBATCH --output=TRIM.%j.out
#################
# a file for errors from the job for each sample
#SBATCH --error=TRIM.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 6:00:00
#################
#quality of service; think of it as job priority
#SBATCH --partition=amilan
#SBATCH --qos=normal
#################
#number of nodes- for RM-shared -N 1
#SBATCH --ntasks-per-node 4
#################
#SBATCH --mem=12G
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=END
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=ldurkee@colostate.edu
#################
#now run normal batch commands
##################

#echo commands to stdout
set -x

##################
#run sbatch in for loop
#AL2_CKDL230002367-1A_HJJGFCCX2_L6_1.fq.gz

#now run normal batch commands within a for loop, run in folder with fastq files to assign correct plate number

#for sample in `ls *1.fq.gz |cut -f1 -d '_'`; do lane = `ls *1.fq.gz |cut -f3 -d '_'`; echo  $sample; sbatch ../1.trimmomatic.sbatch $sample $lane ; done


source ~/.bashrc
module load trimmomatic

##Variables:sampel name and number associated with sample

lane=$2

sample=$1

sample2=$(echo ${sample} | sed 's/1.fq.gz/2.fq.gz/')

id=$(echo ${sample} | cut -f1 -d_)

lane=$(echo ${sample} | cut -f3 -d_)

java -jar $TRIMMOMATIC PE -threads 4 $sample $sample2  \
              "$id"_"$lane"_1.trimmed.fastq "$id"_"$lane"_1un.trimmed.fastq \
              "$id"_"$lane"_2.trimmed.fastq "$id"_"$lane"_2un.trimmed.fastq \
              ILLUMINACLIP:$ADAPTERS/TruSeq3-PE-2.fa:2:30:10 \
              SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:36



