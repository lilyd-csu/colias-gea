## file overview

Follows GATK workflow using GATK4: https://gatk.broadinstitute.org/

5. Creates gVCF file
6. GATK GenomicsDBImport
7. GATK GenotypeGVCFs
8. Removes related individuals using NGSRelate
9. Sets NA values to "missing"
10. Merge and filter to create final VCF
