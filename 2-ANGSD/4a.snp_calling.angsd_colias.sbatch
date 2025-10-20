#!/bin/bash 
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=angsd
#################  
#a file for job output, you can check job progress
#SBATCH --output=angsd.%j.out
#################
# a file for errors from the job
#SBATCH --error=angsd.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 24:00:00
#################
#quality of service; think of it as job priority
#SBATCH -p amilan
#SBATCH --qos=long
#################
#number of nodes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 24
#################
#SBATCH --mem=30G
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

source ~/.bashrc
conda activate angsd

module load jdk
module load python

FASTA="/scratch/alpine/ldurkee@colostate.edu/colias2023/reference_ncbi/reference1/GCF_905220415.1_ilColCroc2.1_genomic.fna"


j=$1 #bam list

outfile=$(echo $j|sed 's/-bam-list.txt/.minInd90.05maf.baq.HiSeq_filt/g')

#NGSTOOLS="/jet/home/cbossu/bridges1-home/programs/ngsTools/angsd/angsd"


#maxdepth=85 (number of individuals)*2 (average coverage)*3 (if you have a lot mapped here, 3x times the ammount, it's likely a repetitive region')
#maxdepth=179*4*3

#This gets you the beagle.gz file for pcangsd
angsd -bam $j -P 24 -ref "$FASTA" -anc "$FASTA" -dosaf 1 \
        -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -baq 1 \
        -maxDepth 2148 -doGlf 2 \
        -minMapQ 20 -minQ 20 -minInd 90 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-6 -doGeno 32 -doPost 1  -minMaf 0.05 -out $outfile


## This gets you the glf.gz file for ngsRelate


