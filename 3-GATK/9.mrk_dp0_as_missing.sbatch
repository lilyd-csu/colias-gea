#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=gvcf
#################
#a file for job output, you can check job progress
#SBATCH --output=gvcf.%j.out
#################
# a file for errors from the job
#SBATCH --error=gvcf.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 8:00:00
#################
#quality of service; think of it as job priority
#SBATCH -p amilan
#################
#number of nodes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 2
#################
#SBATCH --mem=7G
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
##run outside of script for each sample : for vcf in `ls unfilt_genotypes/*vcf.gz|head -1`; do sbatch 7b.mrk_dp0_as_missing.sbatch $vcf; done

source ~/.bashrc
module load bcftools

vcf=$1

##Location of fasta file
FASTA="/scratch/alpine/ldurkee@colostate.edu/colias2023/reference_ncbi/reference1/GCF_905220415.1_ilColCroc2.1_genomic.fna"

#conda activate bcftools

mkdir -p mk_dp0_as_missing/unfilt_genotypes
bcftools +setGT $vcf -- -t q -n . -i 'FMT/DP=0 | (FMT/PL[:0]=0 & FMT/PL[:1]=0 & FMT/PL[:2]=0)' | bcftools +fill-tags - -- -t 'NMISS=N_MISSING' | bcftools view -Oz - > mk_dp0_as_missing/"$vcf"
bcftools index -t mk_dp0_as_missing/"$vcf"
