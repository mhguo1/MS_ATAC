#!/usr/bin/env Rscript
## R script to prepare DICE results for coloc that overlap with MS GWAS (Science 2019 paper)
## Last update: 12/11/2021
## Nikolaos Patsopoulos



## Usage: Rscript /work/Projects/MS-chromatin_fine_mapping/scripts/DICE_prepare_coloc.R -e /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/data_prep/{1}.{2}.complete.txt -o /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/{2}/chr{1}/all.coloc.txt -d /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/{2}/chr{1}

## -- Read options -- ##
if(!require(optparse)){
    install.packages("optparse", repos = "http://cran.us.r-project.org", lib="~/local/R_libs/")
}
library(optparse)



option_list = list(
  make_option(c("-e", "--eqtl"), type="character", default=NULL,
              help="Complete eQTL results (", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="default=NULL",
              help="Dir and prefix of output files", metavar="character"),
  make_option(c("-d", "--dir"), type="character", default="default=NULL",
              help="Dir to output per gene files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## -- ENDOF Read options -- ##


## -- Libraries -- ##
#library(tidyverse)
 

## -- ENDOF Libraries -- ##


## -- Main -- ##
## eqtl 
eqtl <- read.table(opt$eqtl, header=TRUE, stringsAsFactors=FALSE, sep=" ")
#eqtl <- read.table("/work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/data_prep/10.Th17.complete.txt", header=TRUE, stringsAsFactors=FALSE, sep=" ")

## estimate var_beta
eqtl$se_beta <- abs(eqtl$beta / eqtl$z_eqtl)
eqtl$var_beta <- eqtl$se_beta * eqtl$se_beta

## keep needed columns
eqtl2 <- eqtl[ , c(1, 3, 12, 14, 18, 26, 24, 15:17, 20, 21)]

## add sdY
eqtl2$sdY <- "1"

## add type
eqtl2$type <- "quant"

## reformat and rename
eqtl2 <- eqtl2[ , c(1:6, 13:14, 7, 8:12)]
names(eqtl2)[1:6] <- c("gene", "ensembl", "snp", "pos", "beta", "varbeta")

## remove lines with NAs 
eqtl2 <- eqtl2[!is.na(eqtl2$varbeta),]



## output file
write.table(eqtl2, file=opt$out, quote=FALSE, row.names=FALSE, col.names=TRUE)

## output per gene files
# list of genes
genes <- unique(eqtl2$gene)
for( i in 1:length(genes) ){
 print(i)
 # subset summary statistics
 tmp <- eqtl2[  eqtl2$gene==genes[i], ]
 # output
 write.table(tmp, paste0(opt$dir, "/", genes[i], ".coloc.txt"), quote=FALSE, row.names=FALSE)
 
}

## -- ENDOF Main -- ##
