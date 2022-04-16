#!/usr/bin/sh

#Step 1: Munge MS GWAS summary statistics. Summary statistics were obtained from Patsopoulos et al., Science 2019. These were obtained from the International Multiple Sclerosis Genetics Consortium
(https://imsgc.net/)

python /software/ldsc/munge_sumstats.py \
--sumstats ${ms_summary_stats} \
--merge-alleles /software/ldsc_files/w_hm3.snplist \
--out MS \
--nstudy N \
--N 41505 \
--nstudy-min 1 \
--N-cas 14802 \
--N-con 26703

#Step 2: Compute LD scores for ATAC-seq peaks from each cell type
for cell in B CD4 CD8 CLP CMP Ery GMP HSC LMPP MEP MPP Mega Mono NK mDC pDC_ATAC
do
for i in {1..22}
do
python /software/ldsc/make_annot.py \
--bimfile /ldsc_files/1000G_plinkfiles/1000G.mac5eur.${i}.bim \
--bed-file ${cell}_ATAC.peaks.bed.gz \
--annot-file ${cell}_ATAC.${i}.annot.gz

python /software/ldsc/ldsc.py \
--l2 \
--bfile /ldsc_files/1000G_plinkfiles/1000G.mac5eur.${i} \
--ld-wind-cm 1 \
--annot /ldsc_files/cell_type_groups/buenrostro/${cell}_ATAC.${i}.annot.gz \
--out ${cell}_ATAC.${i} \
--print-snps /ldsc_files/hapmap3_snps/hm.${i}.snp
done
done

#Step 3: Partitioned heritability for single annotation
for cell in B CD4 CD8 CLP CMP Ery GMP HSC LMPP MEP MPP Mega Mono NK mDC pDC_ATAC
do

python /cvar/jhlab/mguo/SOFTWARE/ldsc/ldsc.py \
--h2 MS.sumstats.gz \
--overlap-annot \
--ref-ld-chr /ldsc_files/baseline/baseline.,/ldsc_files/cell_type_groups/buenrostro/${cell}_ATAC. \
--w-ld-chr /ldsc_files/weights_hm3_no_hla/weights. \
--frqfile-chr /ldsc_files/1000G_frq/1000G.mac5eur. \
--print-coefficients \
--out ${trait}_buenrostro_ATAC
done
