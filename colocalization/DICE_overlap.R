#!/usr/bin/env Rscript
## R script to find overlap of MS GWAS variants and DICE all.txt results.
## Last update: 12/11/2021
## Nikolaos Patsopoulos



## Usage: Rscript /work/Projects/MS-chromatin_fine_mapping/scripts/DICE_overlap.R -v /work/Projects/MS-chromatin_fine_mapping/results/MS.PICS.ld0.2.ATAC.gencode_v19.PCHiC.motifbreakr.eqtl.cs.bed -e work/Tools/smr/data/DICE/chr10.Th17.all.txt -o /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/data_prep/chr10.Th17


## -- Read options -- ##
if(!require(optparse)){
    install.packages("optparse", repos = "http://cran.us.r-project.org", lib="~/local/R_libs/")
}
library(optparse)



option_list = list(
  make_option(c("-v", "--variants"), type="character", default=NULL,
              help="Text file with PICS+LD info per MS locus", metavar="character"),
  make_option(c("-e", "--eqtl"), type="character", default=NULL,
              help="Formatted eQTL results (all.txt from SMR analyses)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="default=NULL",
              help="Dir and prefix of output files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## -- ENDOF Read options -- ##


## -- Libraries -- ##
library(tidyverse)


## -- ENDOF Libraries -- ##


## -- Main -- ##
## load variants
variants <- read.table(opt$variants, header=TRUE, stringsAsFactors=FALSE, sep="\t")
#variants <- read.table("/work/Projects/MS-chromatin_fine_mapping/results/MS.PICS.ld0.2.ATAC.gencode_v19.PCHiC.motifbreakr.eqtl.cs.bed", header=TRUE, stringsAsFactors=FALSE, sep="\t")

## eqtl
eqtl <- read.table(opt$eqtl, header=TRUE, stringsAsFactors=FALSE, sep=" ")
#eqtl <- read.table("/work/Tools/smr/data/DICE/chr10.Th17.all.txt", header=TRUE, stringsAsFactors=FALSE, sep=" ")
# remove associatiosn with z == 0. these create many issues later on
eqtl <- eqtl[ eqtl$z_eqtl!=0  ,]
## get chromosome
chr <- unique(eqtl$variant_CHR)

## subset
nrow(variants)
variants <- variants[ variants$chr==chr, ]
nrow(variants)


## find overlap
common_v <- variants[ variants$snp_end %in% eqtl$POS_hg19, ]
nrow(common_v)
common_e <- unique(eqtl[ eqtl$POS_hg19 %in% variants$snp_end, 3]) # this is a list of genes
length(common_e)

## extract the list of eqtls
eqtl_c <- eqtl[ eqtl$ensemble_id %in% common_e , ]
nrow(eqtl)
nrow(eqtl_c)

## list of variants that do not overlap
not_overlap_v <- variants[ !(variants$snp_end %in% eqtl$POS_hg19), ]
nrow(not_overlap_v)

## find list of effect-genes
# create light version of the dataframes
v_l <- variants[ , c(1:3, 8)] 
names(v_l)[4] <- "POS_hg19"
e_l <- eqtl[ , c(1, 3, 14 )]
joined <- merge(v_l, e_l)
nrow(joined)
# find the unique effect-gene combinations
joined2 <- joined[ , c(3, 5)] 
nrow(joined2)
joined2 <- unique(joined2)
nrow(joined2)


## output files
# eqtl, this is the primary file
write.table(eqtl_c, file=paste0(opt$out, ".complete.txt"), quote=FALSE, row.names=FALSE)
# common variants
write.table(common_v, file=paste0(opt$out, ".variants.complete.txt"), quote=FALSE, row.names=FALSE)
# non-overlapping variants
write.table(not_overlap_v, file=paste0(opt$out, ".non_overlapping_variants.complete.txt"), quote=FALSE, row.names=FALSE)
# list of effects-genes
write.table(joined2, file=paste0(opt$out, ".effects_genes.txt"), quote=FALSE, row.names=FALSE)



## -- ENDOF Main -- ##
