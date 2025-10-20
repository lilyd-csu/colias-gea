#!/bin/bash 
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=ngsrel
#################  
#a file for job output, you can check job progress
#SBATCH --output=ngsrel.%j.out
#################
# a file for errors from the job
#SBATCH --error=ngsrel.%j.err
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
#SBATCH --ntasks-per-node 10
#################
#SBATCH --mem=20G
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
conda activate binfr

module load jdk
module load python



#BAM1="/scratch/alpine/cbossu@colostate.edu/BUOW/bwa_juv2022"
#COMMS="/scratch/alpine/cbossu@colostate.edu/BUOW/BUOW.j2022.2m.comm_lines.txt"
NGSRELATE="/projects/ldurkee@colostate.edu/ngsRelate"

#My file is called: colias-full-down4.0x.minInd90.05maf.baq.HiSeq_filt_rm.relate.glf.gz
#Path: /scratch/alpine/ldurkee@colostate.edu/colias2023/bwa_mem/angsd_4.0x/

### Then we extract the frequency column from the allele frequency file and 
# remove the header (to make it in the format NgsRelate needs)
zcat colias-full-down4.0x.minInd90.05maf.baq.HiSeq_filt_rm.relate.mafs.gz | cut -f5 |sed 1d > freq

### run NgsRelate
$NGSRELATE/ngsRelate  -g colias-full-down4.0x.minInd90.05maf.baq.HiSeq_filt_rm.relate.glf.gz -n 179 -z colias-full-down4.0x-bam-list.ind -f freq  -O $outfile.newres

