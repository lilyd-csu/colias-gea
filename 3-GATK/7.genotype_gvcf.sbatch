#!/bin/bash 
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=genotype_gvcf
#################  
#a file for job output, you can check job progress
#SBATCH --output=genotype_gvcf.%j.out
#################
# a file for errors from the job
#SBATCH --error=genotype_gvcf.%j.err
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
## 1-1 is test; now run 2-31
#SBATCH --array=2-31
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

source ~/.bashrc
module load gatk


##To prepare the comm_lines.txt file, choose one bam file
#samtools view -H bam_file | sed 's/SN://g; s/LN://g' | awk '/^@SQ/' | awk -f assemble-scaffold-lists.awk > colias.4m.comm_lines.txt

FASTA="/scratch/alpine/ldurkee@colostate.edu/colias2023/reference_ncbi/reference1/GCF_905220415.1_ilColCroc2.1_genomic.fna"
GVCF="/scratch/alpine/ldurkee@colostate.edu/colias2023/bwa_mem/raw_genotypes"
COMMS="/scratch/alpine/ldurkee@colostate.edu/colias2023/bwa_mem/colias.3m.comm_lines.txt"

OUTN=$(printf '%04d' $SLURM_ARRAY_TASK_ID)
OUT=$OUTN

SEG=$(awk -F"\t" -v n=$SLURM_ARRAY_TASK_ID 'NR == n {print $2}' $COMMS)
     
#Step 7: Genotype all GVCF files

mkdir unfilt_genotypes
gatk GenotypeGVCFs \
   -R $FASTA \
   -V gendb://called_genotypes/${OUT}_colias_database \
   -O unfilt_genotypes/${OUT}_colias_genotypes_unfiltered.vcf.gz > unfilt_genotypes/${OUT}.stdout 2>unfilt_genotypes/${OUT}.stderr
