#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=BWA_Dup
#################
#a file for job output, you can check job progress for each sample
#SBATCH --output=BWA-Dup.%j.out
#################
# a file for errors from the job for each sample
#SBATCH --error=BWA-Dup.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 24:00:00
#################
#quality of service; think of it as job priority
#SBATCH --partition=amilan
#SBATCH --qos=normal
#################
#number of nodes- for RM-shared -N 1
#SBATCH --ntasks-per-node 4
#################
#SBATCH --mem=16G
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
#run sbatch in for loop, for L6, L7, and L8 separately
#now run normal batch commands within a for loop, run in folder with fastq files to assign correct plate number
#for sample in `ls *L6_1.trimmed.fastq | cut -f1 -d'_'`; do echo $sample;  sbatch ../2.map_bwa_colias.sbatch $sample L6 ;done

#acompile

module load bwa
module load samtools
module load picard
#create conda environment that has bwa, samtools and picard. Recently I made two conda environments (binfr= bcftools, samtools, bwa; picard= picard)
#source ~/.bashrc
#conda activate binfr

##Variables: Plate number and directory for bamUtil
#PLATE="Plate1"

FASTA="/scratch/alpine/ldurkee@colostate.edu/colias2023/reference_ncbi/reference1/GCF_905220415.1_ilColCroc2.1_genomic.fna"
# reference genome for Colias crocea (clouded yellow)
# European Nucleotide Archive: Accession number PRJEB42949; https://identifiers.org/ena.embl/PRJEB42949

sample=$1
lane=$2

#cd ../
#create index for BUOW.fasta

#run 1-2 below with reference genome
#bwa index  GCF_905220415.1_ilColCroc2.1_genomic.fna
#samtools faidx GCF_905220415.1_ilColCroc2.1_genomic.fna

#cd ../trimmomatic
##Align each sample to genome. Note that genome reference must already be built through bwa
#mkdir ../bwa_mem

#Note, plate is the order the plate was sequenced. 
ID="$lane.$sample"

#MS13_L7_1.trimmed.fastq
#MS13_L7_2.trimmed.fastq

##map paired reads
bwa mem  -t 4 $FASTA "$sample"_"$lane"_1.trimmed.fastq  "$sample"_"$lane"_2.trimmed.fastq > ../bwa_mem/aln_"$sample"_"$lane".sam

cd ../bwa_mem
#########sort, add read group information and index it#########
samtools sort -o aln_"$sample"_"$lane".bam aln_"$sample"_"$lane".sam

##Add read groups
picard AddOrReplaceReadGroups INPUT=aln_"$sample"_"$lane".bam RGID="$ID" RGLB="$lane" RGPL=illumina.HJJGFCCX2 RGPU="$lane"."$sample" \
RGSM="$sample" OUTPUT="$sample"_"$lane"_RG.bam VALIDATION_STRINGENCY=SILENT
samtools index "$sample"_"$lane"_RG.bam

rm aln_"$sample".*

##Removing PCR duplicates from BAM files
#sort by name, not position
samtools sort -n -o "$sample"_"$lane"_namesort.bam "$sample"_"$lane"_RG.bam

#Add mate score tags for samtools markdup to select best reads
samtools fixmate -m "$sample"_"$lane"_namesort.bam "$sample"_"$lane"_fixm.bam
rm "$sample"_"$lane"_namesort.bam
rm "$sample"_"$lane"_RG.bam


#Sort again by position
samtools sort -o "$sample"_"$lane"_fixm.sort.bam  "$sample"_"$lane"_fixm.bam
rm "$sample"_"$lane"_fixm.bam

#Markdups
samtools markdup "$sample"_"$lane"_fixm.sort.bam "$sample"_"$lane"_mkdup.bam
samtools index "$sample"_"$lane"_mkdup.bam
