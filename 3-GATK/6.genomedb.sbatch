#!/bin/bash 
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=genomedb
#################  
#a file for job output, you can check job progress
#SBATCH --output=genomedb.%j.out
#################
# a file for errors from the job
#SBATCH --error=genomedb.%j.err
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
## will do a test with 1-1 - JOB COMPLETED
## now will do 2-31
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
#samtools view -H bam_file | sed 's/SN://g; s/LN://g' | awk '/^@SQ/' | awk -f assemble-scaffold-lists.3Mbp.awk > colias.3m.comm_lines.txt

#change these paths 
FASTA="/scratch/alpine/ldurkee@colostate.edu/colias2023/reference_ncbi/reference1/GCF_905220415.1_ilColCroc2.1_genomic.fna"
COMMS="/scratch/alpine/ldurkee@colostate.edu/colias2023/bwa_mem/colias.3m.comm_lines.txt"

# job1 will be 0001
OUTN=$(printf '%04d' $SLURM_ARRAY_TASK_ID)
OUT=$OUTN

SEG=$(awk -F"\t" -v n=$SLURM_ARRAY_TASK_ID 'NR == n {print $2}' $COMMS)

#6b: Consolidate all GVCF files
#This step consists of consolidating the contents of GVCF files across all of your samples in order to improve scalability and speed the next step, joint genotyping.
# We will use the GenomicsDBImport tool and will need to make a new temporary directory.

mkdir called_genotypes
mkdir called_genotypes/temp

gatk GenomicsDBImport \
   $(printf ' --variant %s ' raw_genotypes/*g.vcf) \
     --genomicsdb-workspace-path called_genotypes/${OUT}_colias_database \
      $SEG > called_genotypes/${OUT}.stdout 2>called_genotypes/${OUT}.stderr
      

