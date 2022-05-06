#!/usr/bin/env Rscript
## R script to perform coloc between MS GWAS (Science 2019) and DICE
## Last update: 12/12/2021
## Nikolaos Patsopoulos



## Usage: Rscript /work/Projects/MS-chromatin_fine_mapping/scripts/coloc_DICE.R -l /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/data_prep/19.nB.effects_genes.txt -v /work/Projects/MS-chromatin_fine_mapping/results/coloc/GWAS/MS/ -e /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/nB/chr19/ -o /work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/nB/chr19/


## -- Read options -- ##
if(!require(optparse)){
    install.packages("optparse", repos = "http://cran.us.r-project.org", lib="~/local/R_libs/")
}
library(optparse)



option_list = list(
  make_option(c("-l", "--list"), type="character", default=NULL,
              help="List of effect-genes", metavar="character"),
  make_option(c("-v", "--variants"), type="character", default=NULL,
              help="Dir where coloc formatted effects exist", metavar="character"),
  make_option(c("-e", "--eqtl"), type="character", default=NULL,
              help="Dir where coloc formatted gene eQTLs exist", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="default=NULL",
              help="Dir and prefix of output files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## -- ENDOF Read options -- ##


## -- Libraries -- ##
#library(tidyverse)
library(stringr)
library(coloc)


## -- ENDOF Libraries -- ##


## -- Main -- ##
## load effect list
l <- read.table(opt$list, header=TRUE, stringsAsFactors=FALSE)
# l <- read.table("/work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/data_prep/19.nB.effects_genes.txt", header=TRUE, stringsAsFactors=FALSE)
# get chr and position in separate columns
l$chr <- str_split_fixed(l$lead_snp, ":", 4)[,1]
l$pos <- str_split_fixed(l$lead_snp, ":", 4)[,2]
# create var for the effect filename
l$effect <- paste0(l$chr, "__", l$pos)

## load paths
variants <- opt$variants
#variants <- "/work/Projects/MS-chromatin_fine_mapping/results/coloc/GWAS/MS/"
eqtl <- opt$eqtl
#eqtl <-"/work/Projects/MS-chromatin_fine_mapping/results/coloc/DICE/nB/chr19/"

## create empty data frame to store results' nrow is equal to the rows of l
df <- data.frame(matrix(ncol = 12, nrow = nrow(l)))

## - loop over effect-gene - ##
for( i in 1:nrow(l) ){
 print(i)
 # effect-gene
 df[ i , 1] <- l$effect[i] 
 # effect
 df[ i , 2] <- l$lead_snp[i] 
 
 # gene
 df[ i , 3] <- l$gene[i] 

 # load effect data
 v_c <- read.table(paste0(variants, l$effect[i], ".coloc.txt"), header=TRUE)   

 # load eqtl data
 e_c <- read.table(paste0(eqtl, l$gene[i], ".coloc.txt"), header=TRUE) 
 # any variants passes FDR of 10%?
 # n of variants in file
 df[ i , 6] <- ifelse(sum(e_c$fdr_eqtl < 0.10)==0, 0, 1)
 
 ## coloc
 # remove rows that will generate NAs in internal coloc functions
 v_c$test <-  pnorm(-abs(v_c$beta/sqrt(v_c$varbeta))) * 2
 v_c <- v_c[!is.na(v_c$test),]
 # n of variants 
 df[ i , 4] <- nrow(v_c)
 e_c$test <-  pnorm(-abs(e_c$beta/sqrt(e_c$varbeta))) * 2
 e_c <- e_c[!is.na(e_c$test),]
 # n of variants 
 df[ i , 5] <- nrow(e_c)

 
 # minumum datasets
 v_c2 <- as.list(v_c[ , c(7, 9, 3, 2)])
 v_c2$type <- v_c$type[1]
 
 names(e_c)[4] <- "position"
 e_c2 <- as.list(e_c[ , c(5:6, 3:4)])
 e_c2$type <- e_c$type[1]
 e_c2$sdY <- e_c$sdY[1]
 
 # run coloc
 res <-  coloc.abf(dataset1=v_c2, dataset2=e_c2)
 # save results
 df[ i , c(7:12)] <- res$summary

}

## add names 
names(df) <- c("effect_ID", "effect", "gene", "n_effect", "n_gene", "Any_snp_pass_fdr_10per", "n_common", "H0", "H1", "H2", "H3", "H4")


## output file
write.table(df, file=paste0(opt$out, "coloc_res.txt"), quote=FALSE, row.names=FALSE)


## -- ENDOF Main -- ##

