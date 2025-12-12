#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=merge
#################
#a file for job output, you can check job progress
#SBATCH --output=merge.%j.out
#################
# a file for errors from the job
#SBATCH --error=merge.%j.err
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

## run this within the unfilt_genotypes folder

set -x

source ~/.bashrc
module load picard
module load jdk
module load python
module load bcftools
module load gatk

##Location of fasta file and
FASTA="/scratch/alpine/ldurkee@colostate.edu/colias2023/reference_ncbi/reference1/GCF_905220415.1_ilColCroc2.1_genomic.fna"


###Merge vcfs of gatk run
picard MergeVcfs \
 $(printf ' -I %s' 0*vcf.gz) \
--OUTPUT colias.4x.merged_gatk.vcf

#conda activate gatk
###Select variants- SNPs
gatk --java-options "-Djava.io.tmpdir=/scratch/alpine/ldurkee@colostate.edu/tmp -Xmx90g" SelectVariants \
     -R $FASTA \
     -V colias.4x.merged_gatk.vcf \
     -select-type SNP \
     --O colias.4x.merged_gatk.SNP.vcf
 

###VCFtools for filtering
vcftools --vcf colias.4x.merged_gatk.SNP.vcf --out colias.4x.merged_gatk.rm.relate.SNP.05maf --min-alleles 2 --max-alleles 2 --max-missing 0.5 --maf 0.05 --max-maf 0.95 --recode --remove related-species.ind

#conda activate gatk

##GATK for hard filtering VQWSR, Auwera et al 2013, with gatk4

gatk --java-options "-Djava.io.tmpdir=/scratch/alpine/ldurkee@colostate.edu/tmp -Xmx90g" VariantFiltration \
-R $FASTA \
-V colias.4x.merged_gatk.rm.relate.SNP.05maf.recode.vcf \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O colias.4x.merged_gatk.rm.relate.SNP.filtered_gatkVQSR2.vcf

#Exclude filtered variants (PASS= pass filters)
gatk --java-options "-Djava.io.tmpdir=/scratch/alpine/ldurkee@colostate.edu/tmp -Xmx90g" SelectVariants \
     -R $FASTA \
     -V colias.4x.merged_gatk.rm.relate.SNP.filtered_gatkVQSR2.vcf \
     --exclude-filtered  true \
     -O colias.4x.merged_gatk.rm.relate.SNP.filtered_gatkVQSR2.PASS.vcf


#conda activate bcftools

##VCFtools for filtering based on missingness, once you remove systematic errors. Play around with this to get the numbers you want and make sense
vcftools --vcf colias.4x.merged_gatk.rm.relate.SNP.filtered_gatkVQSR2.PASS.vcf --out colias.4x.merged_gatk.rm.relate.SNP.filtered_gatkVQSR2.PASS.8miss --max-missing 0.8 --recode 

