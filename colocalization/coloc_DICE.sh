### Prepare DICE eQTL results for coloc with MS GWAS results
## list of effects: we need to make a list of GWAS effects that overlaps OCRs+PCHiC. We can do this per cell type?
## this might be tricky, so we can start by creating a list for the 200 effects or maybe only the ones that the maringal effect was ~ as the conditional. This will result in 200 or less set of files. 
## now for each, we need to identify the region. We can simplify this by adding +-10K (?) around the range of the PICS variantscfor each effect. 
## then, we will have list of regions with the respective top variant.
## for each of these, we can extract the respective regions from the eQTLs. Now, the tricky part is that these might be many genes per cell type. We need to find a rule here...
## this gene list also needs to overlap with the PCHiC per cell type.


  
## we need files with the summary statistics (per chromosome is inf eand can speedup the operations) and the effects (i.e., the top variant per gene)
## DICE provides FDRs. We are going to use any eGene for which the FDR was <10%.
## However, we do not have MAFs. There are few options: i) drop all AT/CG SNPs and take MAF from 1KG, ii) set sdY to 1
## We have already processed files for the SMR analyses and we can reuse most of the files


## Question: should we use the eQTL catalogue harmonized results?? This might make it easier. Perhaps not for this project

## -- setup system
mkdir -p  /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE
cd /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE

screen -S coloc
conda activate smr

## -- step 1. Find overlapping variants
## There is detailed list with any PICS or SNP in LD (r2>0.2). We will use this list to create two set of files: i) complete: these files will capture eGenes that their cis-eQTL genes (FDR<10%) overlap with at least one SNP. ii) refined: these files will subset the eGenes from the complete files so that the GWAS variants will also overlapp OCRs+PCHiCHs for any(?) cell type.  
mkdir data_prep
#parallel -j30 --bar 'Rscript /work/Projects/MS-chromatin_fine_mapping/scripts/DICE_overlap.R -v /work/Projects/MS-chromatin_fine_mapping/results/MS.PICS.ld0.2.ATAC.gencode_v19.PCHiC.motifbreakr.eqtl.cs.bed -e /work/Tools/smr/data/DICE/chr{1}.{2}.all.txt -o /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/data_prep/{1}.{2}' ::: {1..22} ::: Th17 nB TREG_NAIVE TREG_MEM THSTAR TH2 TH1 TFH
parallel -j30 --bar 'Rscript /work/Projects/MS-chromatin_fine_mapping/scripts/DICE_overlap.R -v /work/Projects/MS-chromatin_fine_mapping/results/MS.PICS.ld0.2.ATAC.gencode_v19.PCHiC.motifbreakr.eqtl.cs.bed -e /work/Tools/smr/data/DICE/chr{1}.{2}.all.txt -o /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/data_prep/{1}.{2}' ::: {1..22} ::: NK MONOCYTES M2 CD8_STIM CD8_NAIVE CD4_STIM CD4_NAIVE

## -- step 2. Format file with summary statistics in R and create the effects (SNP-eGene) file
## we need the beta (lnOR) and var of beta, which is te square(^2) of the SE (of lnOR)
## min info: SNP POS beta var_beta N MAF quant/cc
## N is sample size and is used with MAF to estimate sdY for quantitative phenotypes. 
## However, we do not have MAF (for the author provided summary statistics) and it is a bit risky to guestimate the AT/CG SNPs. We are going to fix sdY to 1 and we can perform sensitivities analyses later on with various values.
#for i in Th17 nB TREG_NAIVE TREG_MEM THSTAR TH2 TH1 TFH
for i in NK MONOCYTES M2 CD8_STIM CD8_NAIVE CD4_STIM CD4_NAIVE
do
 mkdir $i
 for j in {1..22}
 do
  mkdir $i/chr$j
 done
done
#parallel -j30 --bar 'Rscript /work/Projects/MS-chromatin_fine_mapping/scripts/DICE_prepare_coloc.R -e /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/data_prep/{1}.{2}.complete.txt -o /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/{2}/chr{1}/all.coloc.txt -d /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/{2}/chr{1}' ::: {1..22} ::: Th17 nB TREG_NAIVE TREG_MEM THSTAR TH2 TH1 TFH
parallel -j30 --bar 'Rscript /work/Projects/MS-chromatin_fine_mapping/scripts/DICE_prepare_coloc.R -e /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/data_prep/{1}.{2}.complete.txt -o /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/{2}/chr{1}/all.coloc.txt -d /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/{2}/chr{1}' ::: {1..22} ::: NK MONOCYTES M2 CD8_STIM CD8_NAIVE CD4_STIM CD4_NAIVE

## -- step 3. Prepare the per effect MS GWAS summary statistics. This can be done only once for all studies

mkdir -p /work/Projects/MS-chromatin_fine_mapping/results/coloc/GWAS/MS

## format file in R
R --no-restore --no-save

# - libraries - #
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(GenomicRanges)
library(AllelicImbalance)
library(dplyr)
library(stringr)

# - load data -#
data <- read.table("/work/Tools/smr/gwas/ms/discovery_metav3.0.meta", header=TRUE)
nrow(data) # 8589719

# - keep only N==15 - #
data <- data[data$N==15, ]
nrow(data) # 6773531

# - remove lines with NAs - #
data <- data[!is.na(data$P),]
nrow(data) # 6773482

# - replace p values == 0 - #
data$P[data$P==0] <- 1e-310

# - beta - #
data$b <- log(data$OR)

# - z score -#
# add z score
data$z <- abs(qnorm( data$P/2, lower.tail=FALSE ))
# add directionality to z scores
data$z[data$b<0] <- -data$z[data$b<0] 

# - se - #
data$se <- data$b/data$z
data <- data[!is.na(data$se),]
nrow(data) # 6773482


# - var - #
data$varbeta <- data$se * data$se
nrow(data) # 6773482

# - sample size - #
data$n <- 41505

#- type - #
data$type <- "cc"

# - add frequencies (US2 WTCCC2 data set, post imputation; same order of alleles) - #
data$freq <- "NA"
# read file with post-imputation frequencies
freq <- read.table(gzfile("/archive/MS/Broad_archive/pca5.newsnp.assoc.dosage.pos.SEcor.txt.gz"), header=TRUE)
# trim
freq <- freq[ , c(1, 5:7)]
# rename alleles, so that we can later do a check
names(freq)[2:3] <- c("A1_fr", "A2_fr")
# merge
all <- merge(data, freq, all.x=TRUE)
nrow(data)
nrow(all)
nrow(all[ all$A1==all$A1_fr, ] ) # all OK!


# - Add rsIDs - #
# create a GRanges object
gr_in <- makeGRangesFromDataFrame(all, keep.extra.columns=TRUE, seqnames.field="CHR", start.field="BP",end.field="BP" )
seqlevelsStyle(gr_in) <- "NCBI"
gr_in2 <- getSnpIdFromLocation(gr_in, SNPlocs.Hsapiens.dbSNP144.GRCh37)
# add SNPs in first granges object. Add and an if else test to compare length of objects(?)
gr_in$SNP_new <- names(gr_in2)
all2 <- as.data.frame(gr_in,  row.names=NULL)
nrow(all2) # 6773153!
# sanity checks. We expect many to have the same ID
nrow(all2[ all2$SNP==all2$SNP_new, ] ) #5619267
nrow(all2[ all2$SNP!=all2$SNP_new, ] ) #1201252

# the SNP_new has some NAs. So, we will populate the SNP name from the SNP col
all2$SNP_new[is.na(all2$SNP_new)] <-  all2$SNP[is.na(all2$SNP_new)]
summary(as.factor(all2$SNP_new)) # no more NAs


# - remove extra rows and format -# 
all2 <- all2[ , c(1, 2, 22, 7:8, 11:12, 14:15, 13, 10, 21, 17) ]
names(all2) <- c("chr", "position", "snp", "A1", "A2", "OR", "beta", "se", "varbeta", "z", "p", "freq", "type")

# - output file for future usage - #
write.table(all2, file="/archive/MS/Broad_archive/discovery_metav3.0.coloc.txt", quote=FALSE, row.names=FALSE)


# - read variants with LD+PICS - #
variants <- read.table("/work/Projects/MS-chromatin_fine_mapping/results/MS.PICS.ld0.2.ATAC.gencode_v19.PCHiC.motifbreakr.eqtl.cs.bed", header=TRUE, sep="\t")

# - find min/max ranges per effect - #
# add the min pos  per effect
variants <- as.data.frame(variants %>%
   group_by(lead_snp) %>%
   mutate(min_pos = min(snp_end)) )
# add the max pos  per effect
variants <- as.data.frame(variants %>%
   group_by(lead_snp) %>%
   mutate(max_pos = max(snp_end)) )

# - create a unique list of the effects - #
effects <- variants[ variants$lead_snp==variants$ld_snp, c(1:2, 6, 8, 59:60)]
effects$chromosome <- as.numeric(str_split_fixed(effects$chr, "chr", 2)[, 2])

# - output per effect files - #
for( i in 1:nrow(effects) ){
 print(i)
 # subset summary statistics
 tmp <- all2[  all2$chr==effects$chromosome[i] & all2$position>(effects$min_pos[i] - 100000) & all2$position<(effects$max_pos[i] + 100000) , ]
 # output
 write.table(tmp, paste0("/work/Projects/MS-chromatin_fine_mapping/results/coloc/GWAS/MS/", effects$chr[i], "__", effects$snp_end[i], ".coloc.txt"), quote=FALSE, row.names=FALSE)
 
}

q()

## the output is one file per effect with the following columns: chr position snp A1 A2 OR beta se varbeta z p freq type

## -- Step 4. run coloc
#parallel -j30 --bar 'Rscript /work/Projects/MS-chromatin_fine_mapping/scripts/coloc_DICE.R -l /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/data_prep/{1}.{2}.effects_genes.txt -v /work/Projects/MS-chromatin_fine_mapping/results/coloc/GWAS/MS/ -e /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/{2}/chr{1}/ -o /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/{2}/chr{1}/' ::: {1..22} ::: Th17 nB TREG_NAIVE TREG_MEM THSTAR TH2 TH1 TFH
parallel -j30 --bar 'Rscript /work/Projects/MS-chromatin_fine_mapping/scripts/coloc_DICE.R -l /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/data_prep/{1}.{2}.effects_genes.txt -v /work/Projects/MS-chromatin_fine_mapping/results/coloc/GWAS/MS/ -e /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/{2}/chr{1}/ -o /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/{2}/chr{1}/' ::: {1..22} ::: NK MONOCYTES M2 CD8_STIM CD8_NAIVE CD4_STIM CD4_NAIVE

## create a joined file
for j in  Th17 nB TREG_NAIVE TREG_MEM THSTAR TH2 TH1 TFH NK MONOCYTES M2 CD8_STIM CD8_NAIVE CD4_STIM CD4_NAIVE
do
 awk 'NR==1 {print "chr", $0}' /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/$j/chr1/coloc_res.txt > /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/$j/coloc_summary.txt
 for i in {1..22}
 do
  awk -v chr=chr$i 'NR>1 {print chr, $0}' /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/$j/chr$i/coloc_res.txt
 done >> /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/$j/coloc_summary.txt
done   

for j in  Th17 nB TREG_NAIVE TREG_MEM THSTAR TH2 TH1 TFH  NK MONOCYTES M2 CD8_STIM CD8_NAIVE CD4_STIM CD4_NAIVE
do
 echo -n $j" "
 awk 'NR>1 && $NF>0.8' /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/$j/coloc_summary.txt | wc -l 
done  
 
## compress
cd /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/
for j in  Th17 nB TREG_NAIVE TREG_MEM THSTAR TH2 TH1 TFH  NK MONOCYTES M2 CD8_STIM CD8_NAIVE CD4_STIM CD4_NAIVE
do
 echo -n 
 tar -czvf $j.tar.gz $j/coloc_summary.txt
done
