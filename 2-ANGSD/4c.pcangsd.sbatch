#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=pcangsd
#################
#a file for job output, you can check job progress
#SBATCH --output=pcangsd.%j.out
#################
# a file for errors from the job
#SBATCH --error=pcangsd.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 2:00:00
#################
#quality of service; think of it as job priority
#SBATCH -p amilan
#################
#number of nodes needs to match # threads
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 24
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

source ~/.bashrc
conda activate pcangsd

module load jdk
module load python

##Traditionally use pcangsd with beagle.gz files
pcangsd -b colias-full-down4.0x-rm.relate.beagle.gz -o colias-full-down4.0x-rm.relate -e 2 -t 24
