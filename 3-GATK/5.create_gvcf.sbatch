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
#SBATCH -t 24:00:00
#################
#quality of service; think of it as job priority
#SBATCH -p amilan
#################
#number of nodes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 24
#################
#SBATCH --mem=90G
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
##run outside of script for each sample
##run for each lane: L5 L6 L7 - tried, times out!
##will try by site
# for sample in `ls x4.0/AL*bam|cut -f 1 -d "_"|head -1`; do sbatch ../5b.create_gvcf.sbatch $sample; done
source ~/.bashrc
module load gatk

sample=$1

##Location of fasta file
FASTA="/scratch/alpine/ldurkee@colostate.edu/colias2023/reference_ncbi/reference1/GCF_905220415.1_ilColCroc2.1_genomic.fna"

#conda activate gatk

###Generate gvcf file


mkdir ../raw_genotypes
#if using gatk as module but need to specify space
gatk --java-options "-Djava.io.tmpdir=/scratch/alpine/ldurkee@colostate.edu/tmp -Xmx90g" HaplotypeCaller \
	-I  ${sample}_L5_mkdup.bam \
	--min-base-quality-score 20 \
	--minimum-mapping-quality 20 \
	--sample-name ${sample} \
	-R $FASTA \
	--emit-ref-confidence GVCF \
	-O ../raw_genotypes/${sample}.raw.snps.indels.g.vcf


