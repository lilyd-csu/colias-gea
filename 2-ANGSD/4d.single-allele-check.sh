#!/bin/bash

#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=HAPLO2
#################
#a file for job output, you can check job progress for each sample
#SBATCH --output=HAPLO2.%j.out
#################
# a file for errors from the job for each sample
#SBATCH --error=HAPLO2.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 12:00:00
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

# load ANGSD
conda init
conda activate angsd

# Run the dual-validation analysis
# -doGlf 2: Probabilistic Beagle format
# -doHaploCall 1: Single-allele sampling
# -doMajorMinor 1: Infers major/minor alleles from Genotype Likelihoods
# -minInd 132: Filters for sites present in at least 80% of individuals

angsd -bam bams.filelist -out colias_validation \
    -GL 1 \
    -doGlf 2 \
    -doHaploCall 1 \
    -doMajorMinor 1 \
    -doMaf 1 \
    -doPost 1 \
    -doCounts 1 \
    -SNP_pval 1e-6 \
    -minMapQ 30 \
    -minQ 20 \
    -minInd 132 \
    -nThreads 12