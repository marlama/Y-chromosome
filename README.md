# Y-chromosome
All steps for Y chromosome control quality and haplogroup inference using Y-Lineage tracker (https://github.com/Shuhua-Group/Y-LineageTracker)

1. We first select only the ISOGG positions used by Y-Lineage Tracker to infer the haplogroups with the command:

    ##hg38
    vcftools --gzvcf yourdata_chrY_hg38.vcf.gz  --positions ISOGG_Positions_chrY.txt --recode --stdout | gzip -c > yourdata_chrY_hg38_ISOGG.vcf.gz 

    ##hg37
    vcftools --gzvcf yourdata_chrY_hg37.vcf.gz  --positions ISOGG_Positions_hg37_chrY.txt --recode --stdout | gzip -c > yourdata_chrY_hg37_ISOGG.vcf.gz 

2. The resultant file goes directly to the control quality scripts YHaplogroups_QC_HG37.sh or YHaplogroups_QC_HG38.sh.

3. After the control quality, we can run the Y-Lineage tracker on the files *._norm_Ref_chr24_UpdateID_Males_SexInfo_hhmiss_LowHet_PopName.vcf.gz

    ##hg38
    LineageTracker classify --vcf yourfile_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_hhmiss_LowHet_PopName.vcf.gz -b 38 --mut-info -a -o yourfile_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_hhmiss_LowHet_PopName

    ##hg37
    LineageTracker classify --vcf yourfile_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_hhmiss_LowHet_PopName.vcf.gz -b 37 --mut-info -a -o yourfile_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_hhmiss_LowHet_PopName
