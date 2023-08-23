#!/bin/bash

#SBATCH --job-name=chrY_MSSNG_db7_2
#SBATCH -N 1 -c 32 --mem 32g --time 1:00:00

module load bcftools plink/2.3dev plink/1.9.beta3a tabix vcftools/0.1.16

#specify DatasetName
pop=MSSNG_db7_2

# Set Dirctory
cd /hpf/largeprojects/tcagstor/users/marlam/ASD/${pop}/chrY/

#Transform to vcf file
#plink2 --bfile /hpf/largeprojects/tcagstor/scratch/marlam/MSSNG/MSSNG_db7/mssng_db7_2_chrY_updateIDs --recode vcf --out ${pop}_chrY_ISOGG

#fit chr
awk '{if($0 !~ /^#/) print "chr"$0; else if (match($0,/(##contig=<ID=)(.*)/,m)) print m[1]"chr"m[2]; else print $0}' ${pop}_chrY_ISOGG.vcf > ${pop}_chrY_ISOGG_chr.vcf

#Index
bgzip ${pop}_chrY_ISOGG_chr.vcf

tabix -p vcf ${pop}_chrY_ISOGG_chr.vcf.gz

#SplitBiallelic Sites
bcftools norm -m -any ${pop}_chrY_ISOGG_chr.vcf.gz --threads 32 -Oz -o ${pop}_chrY_ISOGG_norm.vcf.gz

bcftools index ${pop}_chrY_ISOGG_norm.vcf.gz --threads 32

#Fit Reference Sites
bcftools norm -Oz --check-ref s -f /hpf/largeprojects/tcagstor/users/marlam/ReferenceFiles/hg19.fa --threads 32 -o ${pop}_chrY_ISOGG_norm_Ref.vcf.gz ${pop}_chrY_ISOGG_norm.vcf.gz

bcftools index ${pop}_chrY_ISOGG_norm_Ref.vcf.gz --threads 32

#Transform to Bfile
plink2 --vcf ${pop}_chrY_ISOGG_norm_Ref.vcf.gz --vcf-half-call missing --make-bed --out ${pop}_chrY_ISOGG_norm_Ref

#Update SNP ID
awk '{print  "24""\t"$1"_"$4"\t"$3"\t"$4"\t"$5"\t"$6}' ${pop}_chrY_ISOGG_norm_Ref.bim > ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID.bim

for i in fam bed ; do cp ${pop}_chrY_ISOGG_norm_Ref.${i} ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID.${i} ; done

#Filter Males
plink --bfile ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID --keep ${pop}_Males.txt --make-bed --out ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males

#Filter heterozygous haploid hardcalls and all female chrY 
awk '{print $1"\t"$2"\t"$3"\t"$4"\t""1""\t"$6}' ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males.fam > ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo.fam

for i in bim bed ; do cp ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males.${i} ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo.${i} ; done

plink --bfile ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo --set-hh-missing --make-bed --out ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_miss

#Filter High Heterozygosity
awk '{print "22""\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo.bim > ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_chr22.bim

for i in bed fam ; do cp ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo.${i} ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_chr22.${i} ; done 

plink --bfile ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_chr22 --freqx --out ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_chr22

count=$(wc -l < ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_chr22.fam)

awk -v count="$count" '{if($6/count>=0.1){print $2}}' ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_chr22.frqx > ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_chr22_high_heterozygous_variants.txt

plink --bfile ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_miss --exclude ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_chr22_high_heterozygous_variants.txt --make-bed --out ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_hhmiss_LowHet

#Update Pop Name
awk -v pop="$pop"  '{print $1"\t"$2"\t"pop"\t"$2}' ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_hhmiss_LowHet.fam > ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_hhmiss_LowHet_PopName.txt

plink --bfile ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_hhmiss_LowHet --update-ids ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_hhmiss_LowHet_PopName.txt --make-bed --recode vcf --out ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_hhmiss_LowHet_PopName

#Index
bgzip ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_hhmiss_LowHet_PopName.vcf

tabix -p vcf ${pop}_chrY_ISOGG_norm_Ref_chr24_UpdateID_Males_SexInfo_hhmiss_LowHet_PopName.vcf.gz

exit