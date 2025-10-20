#!/bin/bash

# Define the array of site codes
sites=("BG" "CR" "DL" "FV" "HM" "HT" "KP" "LF" "LL" "MS" "NF" "OC" "SC" "SP" "WY")

# VCF file
vcf_file="colias.4x.merged_gatk.rm.relate.SNP.filtered_gatkVQSR2.PASS.8miss.recode.vcf.gl_impute4.1.vcf.gz"

# Loop through each site code and run the vcftools commands
for site_code in "${sites[@]}"; do
    # Run vcftools for Tajima's D
    vcftools --gzvcf "$vcf_file" --keep "${site_code}.txt" --TajimaD 100000 --out "${site_code}.100kb"
    
    # Run vcftools for window-pi
    vcftools --gzvcf "$vcf_file" --keep "${site_code}.txt" --window-pi 100000 --out "${site_code}.100kb"
done

echo "All processes completed."

