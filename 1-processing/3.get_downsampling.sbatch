#!/bin/bash
#set a job name
#SBATCH --job-name=downsampling
#SBATCH --output=./err-out/downsampling.%A_%a.out
#SBATCH --error=./err-out/downsampling.%A_%a.err
################
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#################
# Note: 4.84G/core or task
#################
#SBATCH --array=3-444 
#one line of the file per line in downsampling_array_full.txt 
#(1-2 is a test for first two lines only)
#444 total lines
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=END
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=ldurkee@colostate.edu
#################

source ~/.bashrc
#conda activate gatk4

module load picard
module load samtools

set -x 

bam=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' downsampling_array_full.txt)
# 20N002019.bam
cov=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $2}' downsampling_array_full.txt)
# x2.0
p=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $3}' downsampling_array_full.txt)
# 0.23

#picard=/projects/mgdesaix@colostate.edu/programs/picard.jar

input=/scratch/alpine/ldurkee@colostate.edu/colias2023/bwa_mem/${bam}
output=/scratch/alpine/ldurkee@colostate.edu/colias2023/bwa_mem/${cov}/${bam}

picard DownsampleSam I=${input} O=${output} P=${p} VALIDATION_STRINGENCY=SILENT
samtools index ${output}

# mkdir -p ./${cov}/depth

# this checks depth
# depth=$(samtools depth -a ${output} | awk '{sum+=$3} END { print "Average = ",sum/NR}')
# echo ${bam} ${depth} >> ${cov}/depth/${cov}.depth.txt
