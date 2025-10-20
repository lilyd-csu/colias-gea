#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=COV
#################
#a file for job output, you can check job progress for each sample
#SBATCH --output=COV.%j.out
#################
# a file for errors from the job for each sample
#SBATCH --error=COV.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 24:00:00
#################
#quality of service; think of it as job priority
#SBATCH --partition=amilan
#SBATCH --qos=normal
#################
#number of nodes or tasks per node
#SBATCH --ntasks-per-node 1
#################
#SBATCH --mem=4G
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

source ~/.bashrc
#conda activate binfr

#WY6_L6_mkdup.bam

mkdir cov
touch cov/covSummary.bedtools.txt
for sample in `ls BG*.bam` ##you don't need to choose one sample, but you do need to choose 10 or less in the ls function
do
        bedtools genomecov -d -ibam $sample   > cov/"$sample"_cov.bedtools.txt
        awk '{if($3<500) {total+=$3; ++lines}} END {print FILENAME," ",total/lines}' cov/"$sample"_cov.bedtools.txt >> cov/covSummary.bedtools.txt
	rm cov/"$sample"_cov.bedtools.txt
done

