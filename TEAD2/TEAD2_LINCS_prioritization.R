## TEAD2 L1000 comparison check calcs 
# Brenna LaBarre
# created May 18, 2020
# last updated Feb 12, 2021


setwd("~/Desktop/MSC_cluedata/TEAD2_prioritized/")

## load libraries
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(qvalue)
library(gplots)


## read in all the files
TEADko <- read.csv("~/Desktop/MSC_cluedata/TEAD2check/TEAD2_crispr.gct.short.csv", header =T, stringsAsFactors = F) # this is actually a knock out list, not crispr
TEADoe <- read.csv("~/Desktop/MSC_cluedata/TEAD2check/TEAD2_over-expression.gct.short.csv", header = T, stringsAsFactors = F)
CD4Tgenes.new <- read.table("~/Dropbox (Partners HealthCare)/MS-chromatin_fine_mapping/results/Fine-mapping/CD4_prioritized_genes.txt", header = F, stringsAsFactors = F)
TEADtft <- read.table("~/Desktop/MSC_cluedata/TEAD2check/TEAD2_TFT.short.txt", header = F, stringsAsFactors = F)
uniqueCD4.new <- read.table("~/Dropbox (Partners HealthCare)/MS-chromatin_fine_mapping/results/Fine-mapping/cd4_unique", header = F, stringsAsFactors = F)
Bgenes.new <- read.table("~/Dropbox (Partners HealthCare)/MS-chromatin_fine_mapping/results/Fine-mapping/B_prioritized_genes.txt", header = F, stringsAsFactors = F)
uniqueB.new <- read.table("~/Dropbox (Partners HealthCare)/MS-chromatin_fine_mapping/results/Fine-mapping/b_unique", header = F, stringsAsFactors = F)
commongenes.new <- read.table("~/Dropbox (Partners HealthCare)/MS-chromatin_fine_mapping/results/Fine-mapping/common_temp", header = F, stringsAsFactors = F)


# Want a distribution of all genes per cell line, rank and in quantiles
dim(TEADko)
dim(TEADoe)

# for the ko data, cell line info are cols 10:18, for OE it's cols 10:17

ko.ordered <- matrix(NA, ncol = 9, nrow = 12328)
for(i in 1:9){
  ko.ordered[order(TEADko[,(i+9)], decreasing = T),i] <- 1:nrow(ko.ordered)
}
ko.ordered <- as.data.frame(ko.ordered)

oe.ordered <- matrix(NA, ncol = 8, nrow = 12328)
for(i in 1:8){
  oe.ordered[order(TEADoe[,(i+9)], decreasing = T),i] <- 1:nrow(oe.ordered)
}
oe.ordered <- as.data.frame(oe.ordered)

#shorten names to just cell line names
colnames(ko.ordered) <- c("A375", "A549", "HA1E", "HCC515", "HEPG2", "HT29", "MCF7", "PC3", "VCAP")
colnames(oe.ordered) <- c("A375", "A549", "HA1E", "HEPG2", "HT29", "MCF7", "PC3", "VCAP")
# add in gene names at end
ko.ordered$gene <- TEADko$pr_gene_symbol
oe.ordered$gene <- TEADoe$pr_gene_symbol


# bin all the genes based on order (rank)
# each matrix has 12328 genes, so thats 1233 genes per bin (last one will have only 1231)
ko.bins <- matrix(nrow = 12328, ncol = 9)
oe.bins <- matrix(nrow = 12328, ncol = 8)

step =  1233
for(i in 1:10){
  if( i <= 9 && i > 1){
    for(j in 1:9){
      ko.bins[ko.ordered[(step*(i-1)):(step*i), j], j] <- i
    }

  }else if(i == 10){
    for(j in 1:9){
      ko.bins[ko.ordered[11095:12328, j], j] <- i
    }
  }else{
    for(j in 1:9){
      ko.bins[ko.ordered[1:step-1, j], j] <- i
    }
  }
}

# same for the OE
step =  1233
for(i in 1:10){
  if( i <= 9 && i > 1){
    for(j in 1:8){
      oe.bins[oe.ordered[(step*(i-1)):(step*i), j], j] <- i
    }

  }else if(i == 10){
    for(j in 1:8){
      oe.bins[oe.ordered[11095:12328, j], j] <- i
    }
  }else{
    for(j in 1:8){
      oe.bins[oe.ordered[1:step-1, j], j] <- i
     }
   }
 }

## okay, now check what bins the genes from our lists fall into


# ## TFT list
dim(ko.bins[which(TEADko[,2] %in% TEADtft[,1]),])
# # of the 1459 things in the TEADtft list, 961 are genes in the ko data 
pdf("TEAD_TFT_quantilespercelline_KO.pdf")
par(mfrow = c(3,3))
for(i in 1:9){
  plot(table(ko.bins[which(TEADko[,2] %in% TEADtft[,1]),i]),
       main = colnames(ko.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off()

pdf("TEAD_TFT_quantilespercelline_OE.pdf")
par(mfrow = c(2,4))
for(i in 1:8){
  plot(table(oe.bins[which(TEADoe[,2] %in% TEADtft[,1]),i]),
       main = colnames(oe.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off()

##CD4 list
dim(ko.bins[which(TEADko[,2] %in% CD4Tgenes.new[,1]),])
# of the 364 things in the CD4Tgenes.new list, 266 are genes in the ko data 

pdf("TEAD_CD4_quantilespercelline_KO.pdf")
par(mfrow = c(3,3))
for(i in 1:9){
  plot(table(ko.bins[which(TEADko[,2] %in% CD4Tgenes.new[,1]),i]), 
       main = colnames(ko.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off() 

pdf("TEAD_CD4_quantilespercelline_OE.pdf")
par(mfrow = c(2, 4))
for(i in 1:8){
  plot(table(oe.bins[which(TEADoe[,2] %in% CD4Tgenes.new[,1]),i]), 
       main = colnames(oe.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off()

## Unique CD4 list 
dim(ko.bins[which(TEADko[,2] %in% uniqueCD4.new[,1]),])
# of the 185 things in the uniqueCD4.new list, 138 are genes in the ko data 

pdf("TEAD_uniqCD4_quantilespercelline_KO.pdf")
par(mfrow = c(3,3))
for(i in 1:9){
  plot(table(ko.bins[which(TEADko[,2] %in% uniqueCD4.new[,1]),i]), 
       main = colnames(ko.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off() 

pdf("TEAD_uniqCD4_quantilespercelline_OE.pdf")
par(mfrow = c(2, 4))
for(i in 1:8){
  plot(table(oe.bins[which(TEADoe[,2] %in% uniqueCD4.new[,1]),i]), 
       main = colnames(oe.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off() 

## Bcell genes
dim(ko.bins[which(TEADko[,2] %in% Bgenes.new[,1]),])
# of the 261 things in the Bgenes.new list, 174 are genes in the ko data 

pdf("TEAD_B_quantilespercelline_KO.pdf")
par(mfrow = c(3,3))
for(i in 1:9){
  plot(table(ko.bins[which(TEADko[,2] %in% Bgenes.new[,1]),i]), 
       main = colnames(ko.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off()

pdf("TEAD_B_quantilespercelline_OE.pdf")
par(mfrow = c(2, 4))
for(i in 1:8){
  plot(table(oe.bins[which(TEADoe[,2] %in% Bgenes.new[,1]),i]), 
       main = colnames(oe.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off()

## unique Bcell genes
dim(ko.bins[which(TEADko[,2] %in% uniqueB.new[,1]),])
# of the 82 things in the uniqueB.new list, 47 are genes in the ko data 

pdf("TEAD_uniqB_quantilespercelline_KO.pdf")
par(mfrow = c(3,3))
for(i in 1:9){
  plot(table(ko.bins[which(TEADko[,2] %in% uniqueB.new[,1]),i]), 
       main = colnames(ko.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off()

pdf("TEAD_uniqB_quantilespercelline_OE.pdf")
par(mfrow = c(2, 4))
for(i in 1:8){
  plot(table(oe.bins[which(TEADoe[,2] %in% uniqueB.new[,1]),i]), 
       main = colnames(oe.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off()

## common genes
dim(ko.bins[which(TEADko[,2] %in% commongenes.new[,1]),])
# of the 178 things in the common genes list, 127 are genes in the ko data 

pdf("TEAD_common_quantilespercelline_KO.pdf")
par(mfrow = c(3,3))
for(i in 1:9){
  plot(table(ko.bins[which(TEADko[,2] %in% commongenes.new[,1]),i]), 
       main = colnames(ko.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off()

pdf("TEAD_common_quantilespercelline_OE.pdf")
par(mfrow = c(2, 4))
for(i in 1:8){
  plot(table(oe.bins[which(TEADoe[,2] %in% commongenes.new[,1]),i]), 
       main = colnames(oe.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off()

## intersection TFT and CD4
intersect(TEADtft[,1], CD4Tgenes.new[,1]) # 38 overlaps
which(TEADko[,2] %in% intersect(TEADtft[,1], CD4Tgenes.new[,1])) #32

pdf("TEAD_CD4inter_quantilespercelline_KO.pdf")
par(mfrow = c(3,3))
for(i in 1:9){
  plot(table(ko.bins[which(TEADko[,2] %in% intersect(TEADtft[,1], CD4Tgenes.new[,1])),i]), 
       main = colnames(ko.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off() 

pdf("TEAD_CD4inter_quantilespercelline_OE.pdf")
par(mfrow = c(2, 4))
for(i in 1:8){
  plot(table(oe.bins[which(TEADoe[,2] %in% intersect(TEADtft[,1], CD4Tgenes.new[,1])),i]), 
       main = colnames(oe.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off()

## intersection TFT and unique CD4
intersect(TEADtft[,1], uniqueCD4.new[,1]) #15
which(TEADko[,2] %in% intersect(TEADtft[,1], uniqueCD4.new[,1])) #12

pdf("TEAD_uniqCD4inter_quantilespercelline_KO.pdf")
par(mfrow = c(3,3))
for(i in 1:9){
  plot(table(ko.bins[which(TEADko[,2] %in% intersect(TEADtft[,1], uniqueCD4.new[,1])),i]), 
       main = colnames(ko.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off()

pdf("TEAD_uniqCD4inter_quantilespercelline_OE.pdf")
par(mfrow = c(2, 4))
for(i in 1:8){
  plot(table(oe.bins[which(TEADoe[,2] %in% intersect(TEADtft[,1], uniqueCD4.new[,1])),i]), 
       main = colnames(oe.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off()

# intersection TFT and Bcell
intersect(TEADtft[,1], Bgenes.new[,1]) # 25 overlaps
which(TEADko[,2] %in% intersect(TEADtft[,1], Bgenes.new[,1])) #21

pdf("TEAD_Binter_quantilespercelline_KO.pdf")
par(mfrow = c(3,3))
for(i in 1:9){
  plot(table(ko.bins[which(TEADko[,2] %in% intersect(TEADtft[,1], Bgenes.new[,1])),i]), 
       main = colnames(ko.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off() 

pdf("TEAD_Binter_quantilespercelline_OE.pdf")
par(mfrow = c(2, 4))
for(i in 1:8){
  plot(table(oe.bins[which(TEADoe[,2] %in% intersect(TEADtft[,1], Bgenes.new[,1])),i]), 
       main = colnames(oe.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off()

# intersection TFT and uniqueB.new
intersect(TEADtft[,1], uniqueB.new[,1]) # 2 overlaps
which(TEADko[,2] %in% intersect(TEADtft[,1], uniqueB.new[,1])) #1

pdf("TEAD_uniqBinter_quantilespercelline_KO.pdf")
par(mfrow = c(3,3))
for(i in 1:9){
  plot(table(ko.bins[which(TEADko[,2] %in% intersect(TEADtft[,1], uniqueB.new[,1])),i]), 
       main = colnames(ko.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off()

pdf("TEAD_uniqBinter_quantilespercelline_OE.pdf")
par(mfrow = c(2, 4))
for(i in 1:8){
  plot(table(oe.bins[which(TEADoe[,2] %in% intersect(TEADtft[,1], uniqueB.new[,1])),i]), 
       main = colnames(oe.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off() 

# intersection TFT and common
intersect(TEADtft[,1], commongenes.new[,1]) # 23 overlaps
which(TEADko[,2] %in% intersect(TEADtft[,1], commongenes.new[,1])) # 20

pdf("TEAD_commoninter_quantilespercelline_KO.pdf")
par(mfrow = c(3,3))
for(i in 1:9){
  plot(table(ko.bins[which(TEADko[,2] %in% intersect(TEADtft[,1], commongenes.new[,1])),i]), 
       main = colnames(ko.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off()

pdf("TEAD_commoninter_quantilespercelline_OE.pdf")
par(mfrow = c(2, 4))
for(i in 1:8){
  plot(table(oe.bins[which(TEADoe[,2] %in% intersect(TEADtft[,1], commongenes.new[,1])),i]), 
       main = colnames(oe.ordered)[i], xlab = "10%-tile", ylab = "count")
}
dev.off()


## create a table with gene|quantile ko|quantile oe for each cell line

# first, make list of unique genes that are in the l1000 set, add the new lists
allTEADgenes.new <-  union(union(union(union(union(TEADtft[,1], uniqueCD4.new[,1]), CD4Tgenes.new[,1]), 
                                       Bgenes.new[,1]), uniqueB.new[,1]), commongenes.new[,1])
# how many are in the L1000 list?
length(which(TEADoe[,2] %in% union(union(union(union(union(TEADtft[,1], uniqueCD4.new[,1]), CD4Tgenes.new[,1]), 
                                               Bgenes.new[,1]), uniqueB.new[,1]), commongenes.new[,1])))
#[1] 1241
allTEAD.L1000genes.new <- which(TEADoe[,2] %in% union(union(union(union(union(TEADtft[,1], uniqueCD4.new[,1]), CD4Tgenes.new[,1]),
                                                                  Bgenes.new[,1]), uniqueB.new[,1]), commongenes.new[,1]))


#initialize a list to store the tables
TEAD.gene.tables.new <- vector(mode = "list", length = 9)

# the oe experiment is missing a cell line so that is accounted for (the 4th cell line in the ko data)
for(i in 1:9){
  if( i <= 3){
    TEAD.gene.tables.new[[i]] <- data.frame(genes = TEADko[allTEAD.L1000genes.new,2],
                                            ko_quant = ko.bins[allTEAD.L1000genes.new,i],
                                            oe_quant = oe.bins[allTEAD.L1000genes.new,i])
    
  }else if(i == 4){
    TEAD.gene.tables.new[[i]] <- data.frame(genes = TEADko[allTEAD.L1000genes.new,2],
                                            ko_quant = ko.bins[allTEAD.L1000genes.new,i])
  }else{
    TEAD.gene.tables.new[[i]] <- data.frame(genes = TEADko[allTEAD.L1000genes.new,2],
                                            ko_quant = ko.bins[allTEAD.L1000genes.new,i],
                                            oe_quant = oe.bins[allTEAD.L1000genes.new,i-1])
  }
  
}

# name each data frame for the cell lines
names(TEAD.gene.tables.new) <- colnames(ko.ordered)[-10]

### save the list 
save(TEAD.gene.tables.new, file = "./genequantbins.prioritized.Rdata")

#Turn these things into individual tables, and for each individual gene list

# Run some statistical tests on the gene lists

## just for TFT genes
# start with an empty matrix to store the results
res <- matrix(nrow=length(TEAD.gene.tables.new), ncol=12)

#loop over the cell lines

for(i in 1:length(TEAD.gene.tables.new)) {
  print(names(TEAD.gene.tables.new)[i])
  res[i, 1] <- names(TEAD.gene.tables.new)[i]
  
  # total number of genes
  res[i, 2] <- length(which(TEADko[,2] %in% TEADtft[,1]))
  
  # number of genes in top/bottom quantiles
  
  n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TEADtft[,1]), 2])[c(1,10)])
  res[i , 3] <- n
  
  #  binomial test
  
  print(binom.test(n, 961, p=0.2))
  res[i, 4] <- binom.test(n, 961, p=0.2)$p.value
  
  #  chi square test
  
  print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TEADtft[,1]) , 2])))
  res[i, 6] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TEADtft[,1]) , 2]))$p.value
  
  # the 4th cell lines doesn't have OE data
  
  if ( i != 4 ) {
    n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TEADtft[,1]), 3])[c(1,10)])
    res[i , 8] <- n
    print(binom.test(n, 961, p=0.2))
    res[i, 9] <- binom.test(n, 961, p=0.2)$p.value
    print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TEADtft[,1]) , 3])) )
    res[i, 11] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TEADtft[,1]) , 3]))$p.value
    
  } else { print("No data")
    res[i, 9] <- "NA"
    res[i, 11] <- "NA"
  }
}

# local FDR

res[, 5] <- qvalue(as.numeric(res[,4]), pi0 = 1)$lfdr
res[, 7] <- qvalue(as.numeric(res[,6]), pi0 = 1)$lfdr
res[, 10] <- qvalue(as.numeric(res[,9]), pi0 = 1)$lfdr
res[, 12] <- qvalue(as.numeric(res[,11]), pi0 = 1)$lfdr

#transform to data.frame

res.tft.df <- as.data.frame(res)

#add names

names(res.tft.df) <- c("Cell_line", "N_genes", "Genes_in_extreme_quantiles_ko", "Binomial_p_ko", "Binomial_lfdr_ko", "Chi_p_ko", "Chi_lfdr_ko",
                       "Genes_in_extreme_quantiles_OE", "Binomial_p_OE", "Binomial_lfdr_OE", "Chi_p_OE", "Chi_lfdr_OE")

#save a text file

write.table(res.tft.df, "TEAD2_L1000_results_TFT.txt", row.names = FALSE, quote = FALSE)



## just for CD4T genes
res <- matrix(nrow=length(TEAD.gene.tables.new), ncol=12)

#loop over the cell lines

for(i in 1:length(TEAD.gene.tables.new)) {
  print(names(TEAD.gene.tables.new)[i])
  res[i, 1] <- names(TEAD.gene.tables.new)[i]
  
  # total number of genes
  res[i, 2] <- length(which(TEADko[,2] %in% CD4Tgenes.new[,1]))
  
  # number of genes in top/low quantiles
  
  n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% CD4Tgenes.new[,1]), 2])[c(1,10)])
  res[i , 3] <- n
  
  #  binomial test
  
  print(binom.test(n, 266, p=0.2))
  res[i, 4] <- binom.test(n, 266, p=0.2)$p.value
  
  #  chi square test
  
  print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% CD4Tgenes.new[,1]) , 2])))
  res[i, 6] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% CD4Tgenes.new[,1]) , 2]))$p.value
  
  # the 4th cell lines doesn't have OE data
  
  if ( i != 4 ) {
    n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% CD4Tgenes.new[,1]), 3])[c(1,10)])
    res[i , 8] <- n
    print(binom.test(n, 266, p=0.2))
    res[i, 9] <- binom.test(n, 266, p=0.2)$p.value
    print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% CD4Tgenes.new[,1]) , 3])) )
    res[i, 11] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% CD4Tgenes.new[,1]) , 3]))$p.value
    
  } else { print("No data")
    res[i, 9] <- "NA"
    res[i, 11] <- "NA"
  }
}

#local FDR

res[, 5] <- qvalue(as.numeric(res[,4]), pi0 = 1)$lfdr
res[, 7] <- qvalue(as.numeric(res[,6]), pi0 = 1)$lfdr
res[, 10] <- qvalue(as.numeric(res[,9]), pi0 = 1)$lfdr
res[, 12] <- qvalue(as.numeric(res[,11]), pi0 = 1)$lfdr

#transform to data.frame

res.CD4T.df <- as.data.frame(res)

#add names

names(res.CD4T.df) <- c("Cell_line", "N_Genes", "Genes_in_extreme_quantiles_ko", "Binomial_p_ko", "Binomial_lfdr_ko", "Chi_p_ko", "Chi_lfdr_ko",
                        "Genes_in_extreme_quantiles_OE", "Binomial_p_OE", "Binomial_lfdr_OE", "Chi_p_OE", "Chi_lfdr_OE")

#save a text file

write.table(res.CD4T.df, "TEAD2_L1000_results_CD4T.txt", row.names = FALSE, quote = FALSE)



## just for unique CD4 genes
res <- matrix(nrow=length(TEAD.gene.tables.new), ncol=12)

#loop over the cell lines

for(i in 1:length(TEAD.gene.tables.new)) {
  print(names(TEAD.gene.tables.new)[i])
  res[i, 1] <- names(TEAD.gene.tables.new)[i]
  
  # total number of genes
  res[i, 2] <- length(which(TEADko[,2] %in% uniqueCD4.new[,1]))
  
  # number of genes in top/low quantiles
  
  n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueCD4.new[,1]), 2])[c(1,10)])
  res[i , 3] <- n
  
  #  binomial test
  
  print(binom.test(n, 139, p=0.2))
  res[i, 4] <- binom.test(n, 139, p=0.2)$p.value
  
  #  chi square test
  
  print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueCD4.new[,1]) , 2])))
  res[i, 6] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueCD4.new[,1]) , 2]))$p.value
  
  # the 4th cell lines doesn't have OE data
  
  if ( i != 4 ) {
    n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueCD4.new[,1]), 3])[c(1,10)])
    res[i , 8] <- n
    print(binom.test(n, 139, p=0.2))
    res[i, 9] <- binom.test(n, 139, p=0.2)$p.value
    print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueCD4.new[,1]) , 3])) )
    res[i, 11] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueCD4.new[,1]) , 3]))$p.value
    
  } else { print("No data")
    res[i, 9] <- "NA"
    res[i, 11] <- "NA"
  }
}

#local FDR


res[, 5] <- qvalue(as.numeric(res[,4]), pi0 = 1)$lfdr
res[, 7] <- qvalue(as.numeric(res[,6]), pi0 = 1)$lfdr
res[, 10] <- qvalue(as.numeric(res[,9]), pi0 = 1)$lfdr
res[, 12] <- qvalue(as.numeric(res[,11]), pi0 = 1)$lfdr

#transform to data.frame

res.uniqCD4.df <- as.data.frame(res)

#add names

names(res.uniqCD4.df) <- c("Cell_line", "N_genes", "Genes_in_extreme_quantiles_ko", "Binomial_p_ko", "Binomial_lfdr_ko", "Chi_p_ko", "Chi_lfdr_ko",
                           "Genes_in_extreme_quantiles_OE", "Binomial_p_OE", "Binomial_lfdr_OE", "Chi_p_OE", "Chi_lfdr_OE")

#save a text file

write.table(res.uniqCD4.df, "TEAD2_L1000_results_uniqCD4.txt", row.names = FALSE, quote = FALSE)


## just for intersection of TFT and unique CD4 genes

tft.uniqCD4.new.intersect <- data.frame(V1 = intersect(TEADtft[,1], uniqueCD4.new[,1]), stringsAsFactors = F)
res <- matrix(nrow=length(TEAD.gene.tables.new), ncol=12)

#loop over the cell lines

for(i in 1:length(TEAD.gene.tables.new)) {
  print(names(TEAD.gene.tables.new)[i])
  res[i, 1] <- names(TEAD.gene.tables.new)[i]
  
  # total number of genes
  res[i, 2] <- length(which(TEADko[,2] %in% tft.uniqCD4.new.intersect[,1]))
  
  # number of genes in top/low quantiles
  
  n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.uniqCD4.new.intersect[,1]), 2])[which(names(table(TEAD.gene.tables.new[[1]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.uniqCD4.new.intersect[,1]),2])) %in% c(1,10))])
  res[i , 3] <- n
  
  #  binomial test
  
  print(binom.test(n, 12, p=0.2))
  res[i, 4] <- binom.test(n, 12, p=0.2)$p.value
  
  #  chi square test
  
  print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.uniqCD4.new.intersect[,1]) , 2])))
  res[i, 6] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.uniqCD4.new.intersect[,1]) , 2]))$p.value
  
  # the 4th cell lines doesn't have OE data
  
  if ( i != 4 ) {
    n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.uniqCD4.new.intersect[,1]), 3])[which(names(table(TEAD.gene.tables.new[[1]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.uniqCD4.new.intersect[,1]),2])) %in% c(1,10))])
    res[i , 8] <- n
    print(binom.test(n, 12, p=0.2))
    res[i, 9] <- binom.test(n, 12, p=0.2)$p.value
    print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.uniqCD4.new.intersect[,1]) , 3])) )
    res[i, 11] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.uniqCD4.new.intersect[,1]) , 3]))$p.value
    
  } else { print("No data")
    res[i, 9] <- "NA"
    res[i, 11] <- "NA"
  }
}

#local FDR


res[, 5] <- qvalue(as.numeric(res[,4]), pi0 = 1)$lfdr
res[, 7] <- qvalue(as.numeric(res[,6]), pi0 = 1)$lfdr
res[, 10] <- qvalue(as.numeric(res[,9]), pi0 = 1)$lfdr
res[, 12] <- qvalue(as.numeric(res[,11]), pi0 = 1)$lfdr

#transform to data.frame

res.TFTuniqCD4.inter.df <- as.data.frame(res)

#add names

names(res.TFTuniqCD4.inter.df) <- c("Cell_line", "N_genes", "Genes_in_extreme_quantiles_ko", "Binomial_p_ko", "Binomial_lfdr_ko", "Chi_p_ko", "Chi_lfdr_ko",
                                    "Genes_in_extreme_quantiles_OE", "Binomial_p_OE", "Binomial_lfdr_OE", "Chi_p_OE", "Chi_lfdr_OE")

#save a text file

write.table(res.TFTuniqCD4.inter.df, "TEAD2_L1000_results_TFTuniqCD4inter.txt", row.names = FALSE, quote = FALSE)



## just for intersection of TFT and CD4 genes
tft.cd4.new.intersect <- data.frame(V1 = intersect(TEADtft[,1], CD4Tgenes.new[,1]), stringsAsFactors = F)

res <- matrix(nrow=length(TEAD.gene.tables.new), ncol=12)

#loop over the cell lines

for(i in 1:length(TEAD.gene.tables.new)) {
  print(names(TEAD.gene.tables.new)[i])
  res[i, 1] <- names(TEAD.gene.tables.new)[i]
  
  # total number of genes
  res[i, 2] <- length(which(TEADko[,2] %in% tft.cd4.new.intersect[,1]))
  
  # number of genes in top/low quantiles
  
  n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.cd4.new.intersect[,1]), 2])[which(names(table(TEAD.gene.tables.new[[1]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.cd4.new.intersect[,1]),2])) %in% c(1,10))])
  res[i , 3] <- n
  
  #  binomial test
  
  print(binom.test(n, 32, p=0.2))
  res[i, 4] <- binom.test(n, 32, p=0.2)$p.value
  
  #  chi square test
  
  print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.cd4.new.intersect[,1]) , 2])))
  res[i, 6] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.cd4.new.intersect[,1]) , 2]))$p.value
  
  # the 4th cell lines doesn't have OE data
  
  if ( i != 4 ) {
    n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.cd4.new.intersect[,1]), 3])[which(names(table(TEAD.gene.tables.new[[1]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.cd4.new.intersect[,1]),2])) %in% c(1,10))])
    res[i , 8] <- n
    print(binom.test(n, 32, p=0.2))
    res[i, 9] <- binom.test(n, 32, p=0.2)$p.value
    print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.cd4.new.intersect[,1]) , 3])) )
    res[i, 11] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.cd4.new.intersect[,1]) , 3]))$p.value
    
  } else { print("No data")
    res[i, 9] <- "NA"
    res[i, 11] <- "NA"
  }
}

#local FDR


res[, 5] <- qvalue(as.numeric(res[,4]), pi0 = 1)$lfdr
res[, 7] <- qvalue(as.numeric(res[,6]), pi0 = 1)$lfdr
res[, 10] <- qvalue(as.numeric(res[,9]), pi0 = 1)$lfdr
res[, 12] <- qvalue(as.numeric(res[,11]), pi0 = 1)$lfdr

#transform to data.frame

res.TFTCD4.new.inter.df <- as.data.frame(res)

#add names

names(res.TFTCD4.new.inter.df) <- c("Cell_line", "N_genes", "Genes_in_extreme_quantiles_ko", "Binomial_p_ko", "Binomial_lfdr_ko", "Chi_p_ko", "Chi_lfdr_ko",
                                    "Genes_in_extreme_quantiles_OE", "Binomial_p_OE", "Binomial_lfdr_OE", "Chi_p_OE", "Chi_lfdr_OE")

#save a text file

write.table(res.TFTCD4.new.inter.df, "./TEAD2_L1000_results_TFTCD4inter.txt", row.names = FALSE, quote = FALSE)



## B genes
res <- matrix(nrow=length(TEAD.gene.tables.new), ncol=12)

#TFT.common.inter <- data.frame(V1 = intersect(TEADtft[,1], commongenes[,1]), stringsAsFactors = F)

#loop over the cell lines

for(i in 1:length(TEAD.gene.tables.new)) {
  print(names(TEAD.gene.tables.new)[i])
  res[i, 1] <- names(TEAD.gene.tables.new)[i]
  
  # total number of genes
  res[i, 2] <- length(which(TEADko[,2] %in% Bgenes.new[,1]))
  
  # number of genes in top/low quantiles
  
  n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% Bgenes.new[,1]), 2])[which(names(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% Bgenes.new[,1]), 2])) %in% c(1,10))])
  res[i , 3] <- n
  
  #  binomial test
  
  print(binom.test(n, 174, p=0.2))
  res[i, 4] <- binom.test(n, 174, p=0.2)$p.value
  
  #  chi square test
  
  print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% Bgenes.new[,1]) , 2])))
  res[i, 6] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% Bgenes.new[,1]) , 2]))$p.value
  
  # the 4th cell lines doesn't have OE data
  
  if ( i != 4 ) {
    n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% Bgenes.new[,1]), 3])[which(names(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% Bgenes.new[,1]), 3])) %in% c(1,10))])
    res[i , 8] <- n
    print(binom.test(n, 174, p=0.2))
    res[i, 9] <- binom.test(n, 174, p=0.2)$p.value
    print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% Bgenes.new[,1]) , 3])) )
    res[i, 11] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% Bgenes.new[,1]) , 3]))$p.value
    
  } else { print("No data")
    res[i, 9] <- "NA"
    res[i, 11] <- "NA"
  }
}

#local FDR


res[, 5] <- qvalue(as.numeric(res[,4]), pi0 = 1)$lfdr
res[, 7] <- qvalue(as.numeric(res[,6]), pi0 = 1)$lfdr
res[, 10] <- qvalue(as.numeric(res[,9]), pi0 = 1)$lfdr
res[, 12] <- qvalue(as.numeric(res[,11]), pi0 = 1)$lfdr

#transform to data.frame

res.Bgenes.new.df <- as.data.frame(res)

#add names

names(res.Bgenes.new.df) <- c("Cell_line", "N_genes", "Genes_in_extreme_quantiles_ko", "Binomial_p_ko", "Binomial_lfdr_ko", "Chi_p_ko", "Chi_lfdr_ko",
                              "Genes_in_extreme_quantiles_OE", "Binomial_p_OE", "Binomial_lfdr_OE", "Chi_p_OE", "Chi_lfdr_OE")

#save a text file

write.table(res.Bgenes.new.df, "./TEAD2_L1000_results_Bgenes.txt", row.names = FALSE, quote = FALSE)


## uniq B
res <- matrix(nrow=length(TEAD.gene.tables.new), ncol=12)

#TFT.common.inter <- data.frame(V1 = intersect(TEADtft[,1], commongenes[,1]), stringsAsFactors = F)

#loop over the cell lines

for(i in 1:length(TEAD.gene.tables.new)) {
  print(names(TEAD.gene.tables.new)[i])
  res[i, 1] <- names(TEAD.gene.tables.new)[i]
  
  # total number of genes
  res[i, 2] <- length(which(TEADko[,2] %in% uniqueB.new[,1]))
  
  # number of genes in top/low quantiles
  
  n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueB.new[,1]), 2])[which(names(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueB.new[,1]), 2])) %in% c(1,10))])
  res[i , 3] <- n
  
  #  binomial test
  
  print(binom.test(n, 47, p=0.2))
  res[i, 4] <- binom.test(n, 47, p=0.2)$p.value
  
  #  chi square test
  
  print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueB.new[,1]) , 2])))
  res[i, 6] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueB.new[,1]) , 2]))$p.value
  
  # the 4th cell lines doesn't have OE data
  
  if ( i != 4 ) {
    n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueB.new[,1]), 3])[which(names(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueB.new[,1]), 3])) %in% c(1,10))])
    res[i , 8] <- n
    print(binom.test(n, 47, p=0.2))
    res[i, 9] <- binom.test(n, 47, p=0.2)$p.value
    print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueB.new[,1]) , 3])) )
    res[i, 11] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueB.new[,1]) , 3]))$p.value
    
  } else { print("No data")
    res[i, 9] <- "NA"
    res[i, 11] <- "NA"
  }
}

#local FDR


res[, 5] <- qvalue(as.numeric(res[,4]), pi0 = 1)$lfdr
res[, 7] <- qvalue(as.numeric(res[,6]), pi0 = 1)$lfdr
res[, 10] <- qvalue(as.numeric(res[,9]), pi0 = 1)$lfdr
res[, 12] <- qvalue(as.numeric(res[,11]), pi0 = 1)$lfdr

#transform to data.frame

res.uniqueB.new.df <- as.data.frame(res)

#add names

names(res.uniqueB.new.df) <- c("Cell_line", "N_genes", "Genes_in_extreme_quantiles_ko", "Binomial_p_ko", "Binomial_lfdr_ko", "Chi_p_ko", "Chi_lfdr_ko",
                               "Genes_in_extreme_quantiles_OE", "Binomial_p_OE", "Binomial_lfdr_OE", "Chi_p_OE", "Chi_lfdr_OE")

#save a text file

write.table(res.uniqueB.new.df, "./TEAD2_L1000_results_uniqueB.txt", row.names = FALSE, quote = FALSE)


## common
res <- matrix(nrow=length(TEAD.gene.tables.new), ncol=12)


#loop over the cell lines

for(i in 1:length(TEAD.gene.tables.new)) {
  print(names(TEAD.gene.tables.new)[i])
  res[i, 1] <- names(TEAD.gene.tables.new)[i]
  
  # total number of genes
  res[i, 2] <- length(which(TEADko[,2] %in% commongenes.new[,1]))
  
  # number of genes in top/low quantiles
  
  n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% commongenes.new[,1]), 2])[which(names(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% commongenes.new[,1]), 2])) %in% c(1,10))])
  res[i , 3] <- n
  
  #  binomial test
  
  print(binom.test(n, 127, p=0.2))
  res[i, 4] <- binom.test(n, 127, p=0.2)$p.value
  
  #  chi square test
  
  print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% commongenes.new[,1]) , 2])))
  res[i, 6] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% commongenes.new[,1]) , 2]))$p.value
  
  # the 4th cell lines doesn't have OE data
  
  if ( i != 4 ) {
    n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% commongenes.new[,1]), 3])[which(names(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% commongenes.new[,1]), 3])) %in% c(1,10))])
    res[i , 8] <- n
    print(binom.test(n, 127, p=0.2))
    res[i, 9] <- binom.test(n, 127, p=0.2)$p.value
    print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% commongenes.new[,1]) , 3])) )
    res[i, 11] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% commongenes.new[,1]) , 3]))$p.value
    
  } else { print("No data")
    res[i, 9] <- "NA"
    res[i, 11] <- "NA"
  }
}

#local FDR

res[, 5] <- qvalue(as.numeric(res[,4]), pi0 = 1)$lfdr
res[, 7] <- qvalue(as.numeric(res[,6]), pi0 = 1)$lfdr
res[, 10] <- qvalue(as.numeric(res[,9]), pi0 = 1)$lfdr
res[, 12] <- qvalue(as.numeric(res[,11]), pi0 = 1)$lfdr

#transform to data.frame

res.commongenes.new.df <- as.data.frame(res)

#add names

names(res.commongenes.new.df) <- c("Cell_line", "N_genes", "Genes_in_extreme_quantiles_ko", "Binomial_p_ko", "Binomial_lfdr_ko", "Chi_p_ko", "Chi_lfdr_ko",
                                   "Genes_in_extreme_quantiles_OE", "Binomial_p_OE", "Binomial_lfdr_OE", "Chi_p_OE", "Chi_lfdr_OE")

#save a text file

write.table(res.commongenes.new.df, "./TEAD2_L1000_results_commongenes.txt", row.names = FALSE, quote = FALSE)


## bgenes inter
res <- matrix(nrow=length(TEAD.gene.tables.new), ncol=12)

TFT.B.inter.new <- data.frame(V1 = intersect(TEADtft[,1], Bgenes.new[,1]), stringsAsFactors = F)

#loop over the cell lines

for(i in 1:length(TEAD.gene.tables.new)) {
  print(names(TEAD.gene.tables.new)[i])
  res[i, 1] <- names(TEAD.gene.tables.new)[i]
  
  # total number of genes
  res[i, 2] <- length(which(TEADko[,2] %in% TFT.B.inter.new[,1]))
  
  # number of genes in top/low quantiles
  
  n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.B.inter.new[,1]), 2])[which(names(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.B.inter.new[,1]), 2])) %in% c(1,10))])
  res[i , 3] <- n
  
  #  binomial test
  
  print(binom.test(n, 21, p=0.2))
  res[i, 4] <- binom.test(n, 21, p=0.2)$p.value
  
  #  chi square test
  
  print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.B.inter.new[,1]) , 2])))
  res[i, 6] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.B.inter.new[,1]) , 2]))$p.value
  
  # the 4th cell lines doesn't have OE data
  
  if ( i != 4 ) {
    n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.B.inter.new[,1]), 3])[which(names(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.B.inter.new[,1]), 3])) %in% c(1,10))])
    res[i , 8] <- n
    print(binom.test(n, 21, p=0.2))
    res[i, 9] <- binom.test(n, 21, p=0.2)$p.value
    print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.B.inter.new[,1]) , 3])) )
    res[i, 11] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.B.inter.new[,1]) , 3]))$p.value
    
  } else { print("No data")
    res[i, 9] <- "NA"
    res[i, 11] <- "NA"
  }
}

#local FDR


res[, 5] <- qvalue(as.numeric(res[,4]), pi0 = 1)$lfdr
res[, 7] <- qvalue(as.numeric(res[,6]), pi0 = 1)$lfdr
res[, 10] <- qvalue(as.numeric(res[,9]), pi0 = 1)$lfdr
res[, 12] <- qvalue(as.numeric(res[,11]), pi0 = 1)$lfdr

#transform to data.frame

res.TFTB.inter.new.df <- as.data.frame(res)

#add names

names(res.TFTB.inter.new.df) <- c("Cell_line", "N_genes", "Genes_in_extreme_quantiles_ko", "Binomial_p_ko", "Binomial_lfdr_ko", "Chi_p_ko", "Chi_lfdr_ko",
                                  "Genes_in_extreme_quantiles_OE", "Binomial_p_OE", "Binomial_lfdr_OE", "Chi_p_OE", "Chi_lfdr_OE")

#save a text file

write.table(res.TFTB.inter.new.df, "./TEAD2_L1000_results_TFTBinter.txt", row.names = FALSE, quote = FALSE)


## uniq b inter
res <- matrix(nrow=length(TEAD.gene.tables.new), ncol=12)

TFT.uniqB.inter.new <- data.frame(V1 = intersect(TEADtft[,1], uniqueB.new[,1]), stringsAsFactors = F)

#loop over the cell lines

for(i in 1:length(TEAD.gene.tables.new)) {
  print(names(TEAD.gene.tables.new)[i])
  res[i, 1] <- names(TEAD.gene.tables.new)[i]
  
  # total number of genes
  res[i, 2] <- length(which(TEADko[,2] %in% TFT.uniqB.inter.new[,1]))
  
  # number of genes in top/low quantiles
  
  n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.uniqB.inter.new[,1]), 2])[which(names(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.uniqB.inter.new[,1]), 2])) %in% c(1,10))])
  res[i , 3] <- n
  
  #  binomial test
  
  print(binom.test(n, 1, p=0.2))
  res[i, 4] <- binom.test(n, 1, p=0.2)$p.value
  
  #  chi square test
  
  if(i != 3){
    print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.uniqB.inter.new[,1]) , 2])))
    res[i, 6] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.uniqB.inter.new[,1]) , 2]))$p.value
  }else{ res[i, 6] <- "NA"} # this is becase for the 3rd cell line both genes fall in same quantile, chi2 wont run
  
  # the 4th cell lines doesn't have OE data
  
  if ( i != 4 ) {
    n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.uniqB.inter.new[,1]), 3])[which(names(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.uniqB.inter.new[,1]), 3])) %in% c(1,10))])
    res[i , 8] <- n
    print(binom.test(n, 1, p=0.2))
    res[i, 9] <- binom.test(n, 1, p=0.2)$p.value
    
    if(i != 3){
      print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.uniqB.inter.new[,1]) , 3])) )
      res[i, 11] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.uniqB.inter.new[,1]) , 3]))$p.value
    }else{res[i, 11] <- "NA"}
    
  } else { print("No data")
    res[i, 9] <- "NA"
    res[i, 11] <- "NA"
  }
}

#local FDR


res[, 5] <- qvalue(as.numeric(res[,4]), pi0 = 1)$lfdr
res[, 7] <- qvalue(as.numeric(res[,6]), pi0 = 1)$lfdr
res[, 10] <- qvalue(as.numeric(res[,9]), pi0 = 1)$lfdr
res[, 12] <- qvalue(as.numeric(res[,11]), pi0 = 1)$lfdr

#transform to data.frame

res.TFTuniqB.inter.new.df <- as.data.frame(res)

#add names

names(res.TFTuniqB.inter.new.df) <- c("Cell_line", "N_genes", "Genes_in_extreme_quantiles_ko", "Binomial_p_ko", "Binomial_lfdr_ko", "Chi_p_ko", "Chi_lfdr_ko",
                                      "Genes_in_extreme_quantiles_OE", "Binomial_p_OE", "Binomial_lfdr_OE", "Chi_p_OE", "Chi_lfdr_OE")

#save a text file

write.table(res.TFTuniqB.inter.new.df, "./TEAD2_L1000_results_TFTuniqBinter.txt", row.names = FALSE, quote = FALSE)





### and common inter
res <- matrix(nrow=length(TEAD.gene.tables.new), ncol=12)

TFT.common.inter.new <- data.frame(V1 = intersect(TEADtft[,1], commongenes.new[,1]), stringsAsFactors = F)

#loop over the cell lines

for(i in 1:length(TEAD.gene.tables.new)) {
  print(names(TEAD.gene.tables.new)[i])
  res[i, 1] <- names(TEAD.gene.tables.new)[i]
  
  # total number of genes
  res[i, 2] <- length(which(TEADko[,2] %in% TFT.common.inter.new[,1]))
  
  # number of genes in top/low quantiles
  
  n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.common.inter.new[,1]), 2])[which(names(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.common.inter.new[,1]), 2])) %in% c(1,10))])
  res[i , 3] <- n
  
  #  binomial test
  
  print(binom.test(n, 20, p=0.2))
  res[i, 4] <- binom.test(n, 20, p=0.2)$p.value
  
  #  chi square test
  
  print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.common.inter.new[,1]) , 2])))
  res[i, 6] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.common.inter.new[,1]) , 2]))$p.value
  
  # the 4th cell lines doesn't have OE data
  
  if ( i != 4 ) {
    n <- sum(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.common.inter.new[,1]), 3])[which(names(table(TEAD.gene.tables.new[[i]][ which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.common.inter.new[,1]), 3])) %in% c(1,10))])
    res[i , 8] <- n
    print(binom.test(n, 20, p=0.2))
    res[i, 9] <- binom.test(n, 20, p=0.2)$p.value
    print(chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.common.inter.new[,1]) , 3])) )
    res[i, 11] <- chisq.test(table(TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.common.inter.new[,1]) , 3]))$p.value
    
  } else { print("No data")
    res[i, 9] <- "NA"
    res[i, 11] <- "NA"
  }
}

#local FDR


res[, 5] <- qvalue(as.numeric(res[,4]), pi0 = 1)$lfdr
res[, 7] <- qvalue(as.numeric(res[,6]), pi0 = 1)$lfdr
res[, 10] <- qvalue(as.numeric(res[,9]), pi0 = 1)$lfdr
res[, 12] <- qvalue(as.numeric(res[,11]), pi0 = 1)$lfdr

#transform to data.frame

res.TFTcommon.inter.new.df <- as.data.frame(res)

#add names

names(res.TFTcommon.inter.new.df) <- c("Cell_line", "N_genes", "Genes_in_extreme_quantiles_ko", "Binomial_p_ko", "Binomial_lfdr_ko", "Chi_p_ko", "Chi_lfdr_ko",
                                       "Genes_in_extreme_quantiles_OE", "Binomial_p_OE", "Binomial_lfdr_OE", "Chi_p_OE", "Chi_lfdr_OE")

#save a text file

write.table(res.TFTcommon.inter.new.df, "./TEAD2_L1000_results_TFTcommoninter.txt", row.names = FALSE, quote = FALSE)



### Okay, now get the gene counts by cell line and gene list
cell.lines.names <- rep(NA, 18)
for(i in names(TEAD.gene.tables.new)){
  cell.lines.names[which(is.na(cell.lines.names))[1]] <- paste(i, "_KO", sep = "")
  cell.lines.names[which(is.na(cell.lines.names))[1]] <- paste(i, "_OE", sep = "")
}


# Starting with TFT
TFT.cellline.quants.new <- matrix(nrow = length(which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% TEADtft[,1])), ncol = 19, data = NA)
TFT.cellline.quants.new[,1] <- as.character(TEAD.gene.tables.new[[1]][which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% TEADtft[,1]),1])
count <- 2
for(i in 1:9){
  if( i != 4){
    TFT.cellline.quants.new[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TEADtft[,1]),2]
    count <- count + 1
    TFT.cellline.quants.new[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TEADtft[,1]),3]
    count <- count + 1
  }else{
    TFT.cellline.quants.new[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TEADtft[,1]),2]
    count <- count + 2
  }
  
}

TFT.cellline.quants.new <- as.data.frame(TFT.cellline.quants.new)
colnames(TFT.cellline.quants.new) <- c("gene",cell.lines.names)

write.table(TFT.cellline.quants.new, "./TFT_genevscell_quantiles.txt", quote = F, row.names = F)


#  now CD4T
CD4T.cellline.quants.new <- matrix(nrow = length(which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% CD4Tgenes.new[,1])), ncol = 19, data = NA)
CD4T.cellline.quants.new[,1] <- as.character(TEAD.gene.tables.new[[1]][which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% CD4Tgenes.new[,1]),1])
count <- 2
for(i in 1:9){
  if( i != 4){
    CD4T.cellline.quants.new[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% CD4Tgenes.new[,1]),2]
    count <- count + 1
    CD4T.cellline.quants.new[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% CD4Tgenes.new[,1]),3]
    count <- count + 1
  }else{
    CD4T.cellline.quants.new[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% CD4Tgenes.new[,1]),2]
    count <- count + 2
  }
  
}

CD4T.cellline.quants.new <- as.data.frame(CD4T.cellline.quants.new)
colnames(CD4T.cellline.quants.new) <- c("gene",cell.lines.names)

write.table(CD4T.cellline.quants.new, "./CD4T_genevscell_quantiles.txt", quote = F, row.names = F)


#  now unique CD4T
uniqCD4.cellline.quants.new <- matrix(nrow = length(which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% uniqueCD4.new[,1])), ncol = 19, data = NA)
uniqCD4.cellline.quants.new[,1] <- as.character(TEAD.gene.tables.new[[1]][which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% uniqueCD4.new[,1]),1])
count <- 2
for(i in 1:9){
  if( i != 4){
    uniqCD4.cellline.quants.new[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueCD4.new[,1]),2]
    count <- count + 1
    uniqCD4.cellline.quants.new[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueCD4.new[,1]),3]
    count <- count + 1
  }else{
    uniqCD4.cellline.quants.new[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueCD4.new[,1]),2]
    count <- count + 2
  }
  
}

uniqCD4.cellline.quants.new <- as.data.frame(uniqCD4.cellline.quants.new)
colnames(uniqCD4.cellline.quants.new) <- c("gene",cell.lines.names)

write.table(uniqCD4.cellline.quants.new, "./uniqCD4_genevscell_quantiles.txt", quote = F, row.names = F)

#  now B genes
Bgenes.new.cellline.quants <- matrix(nrow = length(which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% Bgenes.new[,1])), ncol = 19, data = NA)
Bgenes.new.cellline.quants[,1] <- as.character(TEAD.gene.tables.new[[1]][which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% Bgenes.new[,1]),1])
count <- 2
for(i in 1:9){
  if( i != 4){
    Bgenes.new.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% Bgenes.new[,1]),2]
    count <- count + 1
    Bgenes.new.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% Bgenes.new[,1]),3]
    count <- count + 1
  }else{
    Bgenes.new.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% Bgenes.new[,1]),2]
    count <- count + 2
  }
  
}

Bgenes.new.cellline.quants <- as.data.frame(Bgenes.new.cellline.quants)
colnames(Bgenes.new.cellline.quants) <- c("gene",cell.lines.names)

write.table(Bgenes.new.cellline.quants, "./Bgenes_genevscell_quantiles.txt", quote = F, row.names = F)

# and unique Bgenes.new
uniqueB.new.cellline.quants <- matrix(nrow = length(which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% uniqueB.new[,1])), ncol = 19, data = NA)
uniqueB.new.cellline.quants[,1] <- as.character(TEAD.gene.tables.new[[1]][which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% uniqueB.new[,1]),1])
count <- 2
for(i in 1:9){
  if( i != 4){
    uniqueB.new.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueB.new[,1]),2]
    count <- count + 1
    uniqueB.new.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueB.new[,1]),3]
    count <- count + 1
  }else{
    uniqueB.new.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% uniqueB.new[,1]),2]
    count <- count + 2
  }
  
}

uniqueB.new.cellline.quants <- as.data.frame(uniqueB.new.cellline.quants)
colnames(uniqueB.new.cellline.quants) <- c("gene",cell.lines.names)

write.table(uniqueB.new.cellline.quants, "./uniqueB_genevscell_quantiles.txt", quote = F, row.names = F)

# Finally common genes
commongenes.new.cellline.quants <- matrix(nrow = length(which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% commongenes.new[,1])), ncol = 19, data = NA)
commongenes.new.cellline.quants[,1] <- as.character(TEAD.gene.tables.new[[1]][which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% commongenes.new[,1]),1])
count <- 2
for(i in 1:9){
  if( i != 4){
    commongenes.new.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% commongenes.new[,1]),2]
    count <- count + 1
    commongenes.new.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% commongenes.new[,1]),3]
    count <- count + 1
  }else{
    commongenes.new.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% commongenes.new[,1]),2]
    count <- count + 2
  }
  
}

commongenes.new.cellline.quants <- as.data.frame(commongenes.new.cellline.quants)
colnames(commongenes.new.cellline.quants) <- c("gene",cell.lines.names)

write.table(commongenes.new.cellline.quants, "./commongenes_genevscell_quantiles.txt", quote = F, row.names = F)


# And then we move on to intersections, starting with CD4T
TFT.CD4.inter.cellline.quants.new <- matrix(nrow = length(which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% tft.cd4.new.intersect[,1])), ncol = 19, data = NA)
TFT.CD4.inter.cellline.quants.new[,1] <- as.character(TEAD.gene.tables.new[[1]][which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% tft.cd4.new.intersect[,1]),1])
count <- 2
for(i in 1:9){
  if( i != 4){
    TFT.CD4.inter.cellline.quants.new[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.cd4.new.intersect[,1]),2]
    count <- count + 1
    TFT.CD4.inter.cellline.quants.new[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.cd4.new.intersect[,1]),3]
    count <- count + 1
  }else{
    TFT.CD4.inter.cellline.quants.new[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.cd4.new.intersect[,1]),2]
    count <- count + 2
  }
  
}

TFT.CD4.inter.cellline.quants.new <- as.data.frame(TFT.CD4.inter.cellline.quants.new)
colnames(TFT.CD4.inter.cellline.quants.new) <- c("gene",cell.lines.names)

write.table(TFT.CD4.inter.cellline.quants.new, "./TFT.CD4.inter_genevscell_quantiles.txt", quote = F, row.names = F)

## and then uniqCD4
tft.uniqCD4.new.intersect.cellline.quants <- matrix(nrow = length(which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% tft.uniqCD4.new.intersect[,1])), ncol = 19, data = NA)
tft.uniqCD4.new.intersect.cellline.quants[,1] <- as.character(TEAD.gene.tables.new[[1]][which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% tft.uniqCD4.new.intersect[,1]),1])
count <- 2
for(i in 1:9){
  if( i != 4){
    tft.uniqCD4.new.intersect.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.uniqCD4.new.intersect[,1]),2]
    count <- count + 1
    tft.uniqCD4.new.intersect.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.uniqCD4.new.intersect[,1]),3]
    count <- count + 1
  }else{
    tft.uniqCD4.new.intersect.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% tft.uniqCD4.new.intersect[,1]),2]
    count <- count + 2
  }
  
}

tft.uniqCD4.new.intersect.cellline.quants <- as.data.frame(tft.uniqCD4.new.intersect.cellline.quants)
colnames(tft.uniqCD4.new.intersect.cellline.quants) <- c("gene",cell.lines.names)

write.table(tft.uniqCD4.new.intersect.cellline.quants, "./TFT.uniqCD4.inter_genevscell_quantiles.txt", quote = F, row.names = F)


## then b cells
TFT.B.inter.new.cellline.quants <- matrix(nrow = length(which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% TFT.B.inter.new[,1])), ncol = 19, data = NA)
TFT.B.inter.new.cellline.quants[,1] <- as.character(TEAD.gene.tables.new[[1]][which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% TFT.B.inter.new[,1]),1])
count <- 2
for(i in 1:9){
  if( i != 4){
    TFT.B.inter.new.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.B.inter.new[,1]),2]
    count <- count + 1
    TFT.B.inter.new.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.B.inter.new[,1]),3]
    count <- count + 1
  }else{
    TFT.B.inter.new.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.B.inter.new[,1]),2]
    count <- count + 2
  }
  
}

TFT.B.inter.new.cellline.quants <- as.data.frame(TFT.B.inter.new.cellline.quants)
colnames(TFT.B.inter.new.cellline.quants) <- c("gene",cell.lines.names)

write.table(TFT.B.inter.new.cellline.quants, "./TFT.B.inter_genevscell_quantiles.txt", quote = F, row.names = F)


## uniq B 
TFT.uniqB.inter.new.cellline.quants <- matrix(nrow = length(which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% TFT.uniqB.inter.new[,1])), ncol = 19, data = NA)
TFT.uniqB.inter.new.cellline.quants[,1] <- as.character(TEAD.gene.tables.new[[1]][which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% TFT.uniqB.inter.new[,1]),1])
count <- 2
for(i in 1:9){
  if( i != 4){
    TFT.uniqB.inter.new.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.uniqB.inter.new[,1]),2]
    count <- count + 1
    TFT.uniqB.inter.new.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.uniqB.inter.new[,1]),3]
    count <- count + 1
  }else{
    TFT.uniqB.inter.new.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.uniqB.inter.new[,1]),2]
    count <- count + 2
  }
  
}

TFT.uniqB.inter.new.cellline.quants <- as.data.frame(TFT.uniqB.inter.new.cellline.quants)
colnames(TFT.uniqB.inter.new.cellline.quants) <- c("gene",cell.lines.names)

write.table(TFT.uniqB.inter.new.cellline.quants, "./TFT.uniqB.inter_genevscell_quantiles.txt", quote = F, row.names = F)


## and then the common genes again 
TFT.commongenes.new.inter.cellline.quants <- matrix(nrow = length(which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% TFT.common.inter.new[,1])), ncol = 19, data = NA)
TFT.commongenes.new.inter.cellline.quants[,1] <- as.character(TEAD.gene.tables.new[[1]][which(as.character(TEAD.gene.tables.new[[1]][,1]) %in% TFT.common.inter.new[,1]),1])
count <- 2
for(i in 1:9){
  if( i != 4){
    TFT.commongenes.new.inter.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.common.inter.new[,1]),2]
    count <- count + 1
    TFT.commongenes.new.inter.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.common.inter.new[,1]),3]
    count <- count + 1
  }else{
    TFT.commongenes.new.inter.cellline.quants[,count] <- TEAD.gene.tables.new[[i]][which(as.character(TEAD.gene.tables.new[[i]][,1]) %in% TFT.common.inter.new[,1]),2]
    count <- count + 2
  }
  
}

TFT.commongenes.new.inter.cellline.quants <- as.data.frame(TFT.commongenes.new.inter.cellline.quants)
colnames(TFT.commongenes.new.inter.cellline.quants) <- c("gene",cell.lines.names)

write.table(TFT.commongenes.new.inter.cellline.quants, "./TFT.commongenes.inter_genevscell_quantiles.txt", quote = F, row.names = F)


### okay, now to create some figures 

# heatmaps of all the quantile counts per cell line/experiment

TEAD2_forheatmaps.new <- list()
count <- 1
ko.cols <- c(2,4,6,8,10,12,14,16,18)
oe.cols <- c(3,5,7,11,13,15,17,19)
for( i in 1:length(ls()[intersect(grep("new", ls()), grep("cellline", ls()))])){
  ref.mat <- get(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][i])
  
  TEAD2_forheatmaps.new[[count]] <- apply(ref.mat[,ko.cols], 2, table)
  count <- count + 1
  
  TEAD2_forheatmaps.new[[count]] <- apply(ref.mat[,oe.cols], 2, table)
  count <- count + 1
}

names(TEAD2_forheatmaps.new) <- c(paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][1], "ko", sep = "_"), paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][1], "oe", sep = "_"),
                                  paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][2], "ko", sep = "_"), paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][2], "oe", sep = "_"),
                                  paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][3], "ko", sep = "_"), paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][3], "oe", sep = "_"),
                                  paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][4], "ko", sep = "_"), paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][4], "oe", sep = "_"),
                                  paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][5], "ko", sep = "_"), paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][5], "oe", sep = "_"),
                                  paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][6], "ko", sep = "_"), paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][6], "oe", sep = "_"),
                                  paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][7], "ko", sep = "_"), paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][7], "oe", sep = "_"),
                                  paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][8], "ko", sep = "_"), paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][8], "oe", sep = "_"),
                                  paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][9], "ko", sep = "_"), paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][9], "oe", sep = "_"),
                                  paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][10], "ko", sep = "_"), paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][10], "oe", sep = "_"),
                                  paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][11], "ko", sep = "_"), paste(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][11], "oe", sep = "_"))

# apply(table) isn't working for some of the data sets, and for the rest the 10 grouping is 
# being alphabetized to row 2, so I have to go through and fix those

#1:6 okay, 7:10 bad, 11:12 good, 13:18 bad, 19:22 good

for(i in c(1,2,3,4,5,6,11,12,19,20,21,22)){
  TEAD2_forheatmaps.new[[i]] <- TEAD2_forheatmaps.new[[i]][c(1,3:10,2),]
}


## The problem with the apply is that not all of them have the same # of quantiles
class(TEAD2_forheatmaps.new[[1]])

for(i in 1:22){
  print(class(TEAD2_forheatmaps.new[[i]]))
}

# 7:10,13:18 are not currently matrices, odds have 9 cols, evens 8 (because OE is missing a cell line) 

for(i in c(7,9,13,15,17)){
  TEAD2_forheatmaps.new[[i]] <- matrix(nrow = 10, ncol = 9, data = NA)
  colnames(TEAD2_forheatmaps.new[[i]]) <- colnames(TEAD2_forheatmaps.new[[1]])
  rownames(TEAD2_forheatmaps.new[[i]]) <- 1:10
}

for(i in c(8,10, 14, 16, 18)){
  TEAD2_forheatmaps.new[[i]] <- matrix(nrow = 10, ncol = 8, data = NA)
  colnames(TEAD2_forheatmaps.new[[i]]) <- colnames(TEAD2_forheatmaps.new[[2]])
  rownames(TEAD2_forheatmaps.new[[i]]) <- 1:10
}

# Now fill in those missing ones (7:10,13:18)
# which are b.inter ko oe, cd4inter ko oe, common inter ko oe, unibinter ko oe, uniqcd4inter ko oe
# basically run the same as above, but manually fill in each row

ls()[intersect(grep("new", ls()), grep("cellline", ls()))][4] # binter 
ref.mat <- get(ls()[intersect(grep("new", ls()), grep("cellline", ls()))][4])
apply(ref.mat[,ko.cols], 2, table)
apply(ref.mat[,oe.cols], 2, table)

TEAD2_forheatmaps.new[[7]][,1] <- c(1,3,3,4,4,NA,NA,2,3,1)
TEAD2_forheatmaps.new[[7]][,2] <- c(3,1,3,3,NA,3,3,NA,2,3)
TEAD2_forheatmaps.new[[7]][,3] <- c(3,1,1,1,3,3,4,1,3,1)
TEAD2_forheatmaps.new[[7]][,4] <- c(2,1,2,5,3,3,1,1,NA,3)
TEAD2_forheatmaps.new[[7]][,5] <- c(4,NA,2,3,1,2,2,4,2,1)
TEAD2_forheatmaps.new[[7]][,6] <- c(NA,1,1,2,4,NA,3,4,3,3)
TEAD2_forheatmaps.new[[7]][,7] <- c(3,2,2,1,3,3,2,3,2,NA)
TEAD2_forheatmaps.new[[7]][,8] <- c(4,3,2,2,NA,3,1,3,2,1)
TEAD2_forheatmaps.new[[7]][,9] <- c(2,NA,6,2,2,1,5,2,1,NA)

TEAD2_forheatmaps.new[[8]][,1] <- c(2,1,2,3,2,4,1,2,3,1)
TEAD2_forheatmaps.new[[8]][,2] <- c(2,4,7,1,3,2,1,NA,1,NA)
TEAD2_forheatmaps.new[[8]][,3] <- c(4,2,NA,4,1,2,2,1,2,3)
TEAD2_forheatmaps.new[[8]][,4] <- c(3,2,6,2,1,4,1,NA,1,1)
TEAD2_forheatmaps.new[[8]][,5] <- c(2,2,5,2,1,2,1,1,2,3)
TEAD2_forheatmaps.new[[8]][,6] <- c(2,2,3,2,6,2,1,NA,2,1)
TEAD2_forheatmaps.new[[8]][,7] <- c(1,1,1,2,1,3,2,5,1,4)
TEAD2_forheatmaps.new[[8]][,8] <- c(2,3,2,3,3,1,1,1,3,2)

ls()[intersect(grep("new", ls()), grep("cellline", ls()))][5] # cd4inter
TEAD2_forheatmaps.new[[9]][,1] <- c(1,5,4,3,6,4,NA,3,4,2)
TEAD2_forheatmaps.new[[9]][,2] <- c(3,2,5,5,NA,3,3,1,5,5)
TEAD2_forheatmaps.new[[9]][,3] <- c(4,1,1,3,7,3,5,2,5,1)
TEAD2_forheatmaps.new[[9]][,4] <- c(2,3,2,8,4,4,2,3,NA,4)
TEAD2_forheatmaps.new[[9]][,5] <- c(4,2,5,3,NA,3,4,5,4,2)
TEAD2_forheatmaps.new[[9]][,6] <- c(1,3,2,3,6,1,5,3,4,4)
TEAD2_forheatmaps.new[[9]][,7] <- c(2,3,3,2,3,4,6,5,3,1)
TEAD2_forheatmaps.new[[9]][,8] <- c(4,4,5,4,1,2,1,4,2,5)
TEAD2_forheatmaps.new[[9]][,9] <- c(3,1,6,4,2,4,6,2,2,2)

TEAD2_forheatmaps.new[[10]][,1] <- c(3,2,3,4,5,5,3,2,6,1)
TEAD2_forheatmaps.new[[10]][,2] <- c(2,4,9,3,4,4,2,1,3,NA)
TEAD2_forheatmaps.new[[10]][,3] <- c(4,5,1,5,1,2,4,3,3,4)
TEAD2_forheatmaps.new[[10]][,4] <- c(4,4,5,4,2,5,4,1,1,2)
TEAD2_forheatmaps.new[[10]][,5] <- c(5,3,5,3,3,4,1,4,1,3)
TEAD2_forheatmaps.new[[10]][,6] <- c(3,4,3,4,6,3,2,1,3,3)
TEAD2_forheatmaps.new[[10]][,7] <- c(1,2,2,5,2,6,2,7,NA,5)
TEAD2_forheatmaps.new[[10]][,8] <- c(2,3,5,5,4,2,2,1,5,3)



ls()[intersect(grep("new", ls()), grep("cellline", ls()))][7] # common inter
TEAD2_forheatmaps.new[[13]][,1] <- c(1,3,3,3,4,NA,NA,2,3,1)
TEAD2_forheatmaps.new[[13]][,2] <- c(2,1,3,3,NA,3,3,NA,2,3)
TEAD2_forheatmaps.new[[13]][,3] <- c(3,1,1,1,3,3,3,1,3,1)
TEAD2_forheatmaps.new[[13]][,4] <- c(2,1,1,5,3,3,1,1,NA,3)
TEAD2_forheatmaps.new[[13]][,5] <- c(4,NA,2,3,NA,2,2,4,2,1)
TEAD2_forheatmaps.new[[13]][,6] <- c(NA,1,1,2,4,NA,3,3,3,3)
TEAD2_forheatmaps.new[[13]][,7] <- c(2,2,2,1,3,3,2,3,2,NA)
TEAD2_forheatmaps.new[[13]][,8] <- c(4,3,2,2,NA,2,1,3,2,1)
TEAD2_forheatmaps.new[[13]][,9] <- c(2,NA,6,2,2,1,4,2,1,NA)

TEAD2_forheatmaps.new[[14]][,1] <- c(2,1,2,3,2,3,1,2,3,1)
TEAD2_forheatmaps.new[[14]][,2] <- c(2,3,7,1,3,2,1,NA,1,NA)
TEAD2_forheatmaps.new[[14]][,3] <- c(4,2,NA,4,NA,2,2,1,2,3)
TEAD2_forheatmaps.new[[14]][,4] <- c(3,2,5,2,1,4,1,NA,1,1)
TEAD2_forheatmaps.new[[14]][,5] <- c(2,2,5,2,1,2,1,1,1,3)
TEAD2_forheatmaps.new[[14]][,6] <- c(2,2,3,2,5,2,1,NA,2,1)
TEAD2_forheatmaps.new[[14]][,7] <- c(1,1,1,2,1,3,2,5,NA,4)
TEAD2_forheatmaps.new[[14]][,8] <- c(2,3,1,3,3,1,1,1,3,2)

ls()[intersect(grep("new", ls()), grep("cellline", ls()))][9] #uniqCD4 inter
TEAD2_forheatmaps.new[[17]][,1] <- c(NA,2,1,NA,2,4,NA,1,1,1)
TEAD2_forheatmaps.new[[17]][,2] <- c(1,1,2,2,NA,NA,NA,1,3,2)
TEAD2_forheatmaps.new[[17]][,3] <- c(1,NA,NA,2,4,NA,2,1,2,NA)
TEAD2_forheatmaps.new[[17]][,4] <- c(NA,2,1,3,1,1,1,2,NA,1)
TEAD2_forheatmaps.new[[17]][,5] <- c(NA,2,3,NA,NA,1,2,1,2,1)
TEAD2_forheatmaps.new[[17]][,6] <- c(1,2,1,1,2,1,2,NA,1,1)
TEAD2_forheatmaps.new[[17]][,7] <- c(NA,1,1,1,NA,1,4,2,1,1)
TEAD2_forheatmaps.new[[17]][,8] <- c(NA,1,3,2,1,NA,NA,1,NA,4)
TEAD2_forheatmaps.new[[17]][,9] <- c(1,1,NA,2,NA,3,2,NA,1,2)

TEAD2_forheatmaps.new[[18]][,1] <- c(1,1,1,1,1,2,2,NA,3,NA)
TEAD2_forheatmaps.new[[18]][,2] <- c(NA,1,2,2,1,2,1,1,2,NA)
TEAD2_forheatmaps.new[[18]][,3] <- c(NA,3,1,1,1,NA,2,2,1,1)
TEAD2_forheatmaps.new[[18]][,4] <- c(1,2,NA,2,1,1,3,1,NA,1)
TEAD2_forheatmaps.new[[18]][,5] <- c(3,1,NA,1,2,2,NA,3,NA,NA)
TEAD2_forheatmaps.new[[18]][,6] <- c(1,2,NA,2,1,1,1,1,1,2)
TEAD2_forheatmaps.new[[18]][,7] <- c(NA,1,1,3,1,3,NA,2,NA,1)
TEAD2_forheatmaps.new[[18]][,8] <- c(NA,NA,4,2,1,1,1,NA,2,1)

ls()[intersect(grep("new", ls()), grep("cellline", ls()))][8] #uniqb inter
TEAD2_forheatmaps.new[[15]][,1] <- c(NA,NA,NA,1,NA,NA,NA,NA,NA,NA)
TEAD2_forheatmaps.new[[15]][,2] <- c(1,NA,NA,NA,NA,NA,NA,NA,NA,NA)
TEAD2_forheatmaps.new[[15]][,3] <- c(NA,NA,NA,NA,NA,NA,1,NA,NA,NA)
TEAD2_forheatmaps.new[[15]][,4] <- c(NA,NA,1,NA,NA,NA,NA,NA,NA,NA)
TEAD2_forheatmaps.new[[15]][,5] <- c(NA,NA,NA,NA,1,NA,NA,NA,NA,NA)
TEAD2_forheatmaps.new[[15]][,6] <- c(NA,NA,NA,NA,NA,NA,NA,1,NA,NA)
TEAD2_forheatmaps.new[[15]][,7] <- c(1,NA,NA,NA,NA,NA,NA,NA,NA,NA)
TEAD2_forheatmaps.new[[15]][,8] <- c(NA,NA,NA,NA,NA,1,NA,NA,NA,NA)
TEAD2_forheatmaps.new[[15]][,9] <- c(NA,NA,NA,NA,NA,NA,1,NA,NA,NA)

TEAD2_forheatmaps.new[[16]][,1] <- c(NA,NA,NA,NA,NA,1,NA,NA,NA,NA)
TEAD2_forheatmaps.new[[16]][,2] <- c(NA,1,NA,NA,NA,NA,NA,NA,NA,NA)
TEAD2_forheatmaps.new[[16]][,3] <- c(NA,NA,NA,NA,1,NA,NA,NA,NA,NA)
TEAD2_forheatmaps.new[[16]][,4] <- c(NA,NA,1,NA,NA,NA,NA,NA,NA,NA)
TEAD2_forheatmaps.new[[16]][,5] <- c(NA,NA,NA,NA,NA,NA,NA,NA,1,NA)
TEAD2_forheatmaps.new[[16]][,6] <- c(NA,NA,NA,NA,1,NA,NA,NA,NA,NA)
TEAD2_forheatmaps.new[[16]][,7] <- c(NA,NA,NA,NA,NA,NA,NA,NA,1,NA)
TEAD2_forheatmaps.new[[16]][,8] <- c(NA,NA,1,NA,NA,NA,NA,NA,NA,NA)

# now to make the heatmaps

tead.cellgene.order <- c("B cells", "CD4+ T cells", "Common", "TFT and B intersection",
                         "TFT and CD4+ T intersection", "TFT", "TFT and common intersection",
                         "TFT and unique to B cell intersection", "TFT and unique to CD4+ T cell intersection",
                         "Unique to CD4+ T cells", "Unique to B cells")
count <- 1
for(i in 1:11){
  if(i != 8){ print(i)
    png(paste("./", ls()[intersect(grep("new", ls()), grep("cellline", ls()))][i], "_heatmap.png", sep = ""), width = 800, height = 600)
    heatmap.2(cbind(TEAD2_forheatmaps.new[[count]], TEAD2_forheatmaps.new[[(count + 1)]]), 
              Rowv = NA,  ylab = "10%-tile",
              main = paste("Genes per Quantile for", tead.cellgene.order[i], sep = " "), col = brewer.pal(10,"RdBu"),
              trace = "none", key = T, density.info = "none", dendrogram = "column", scale = "none")
    dev.off()
    count <- count + 2}
  else{
    count <- count + 2
  }
}



# tft uniq b inter isnt working (BECAUSE THERE IS ONLY ONE GENE)
heatmap.2(cbind(TEAD2_forheatmaps.new[[15]], TEAD2_forheatmaps.new[[16]]), 
          Rowv = FALSE,  ylab = "10%-tile",
          main = paste("Genes per Quantile for", tead.cellgene.order[8], sep = " "), col = brewer.pal(10,"RdBu"),
          trace = "none", key = T, density.info = "none", dendrogram = "none", scale = "none")

# Try heatmaps with the ranks of the genes, instead of the counts per quantile

# need to organize the data a little too
# interested in ko.ordered and oe.ordered data frames

# first, tft

# combine the data
colnames(TEADtft) <- "pr_gene_symbol"
head(merge(TEADtft, TEADko))
merge.tft.ko <- merge(TEADtft, TEADko[,c(2,10:18)])
merge.tft.ko.oe <- merge(merge.tft.ko, TEADoe[,c(2,10:17)], by = "pr_gene_symbol")

# initialize data.frame
TFT.generank.forheatmap.new <- data.frame(gene = merge.tft.ko.oe[,1], stringsAsFactors = F)


for(i in 1:17){
  print(i)
  if(i <= 9){
    TFT.generank.forheatmap.new[,(i+1)] <- order(merge.tft.ko.oe[,(i+1)])
    colnames(TFT.generank.forheatmap.new)[(i+1)] <- paste(colnames(ko.ordered)[i], "_ko", sep = "")
    
  }else{
    TFT.generank.forheatmap.new[,(i+1)] <- order(merge.tft.ko.oe[,(i+1)])
    colnames(TFT.generank.forheatmap.new)[(i+1)] <- paste(colnames(oe.ordered)[(i-9)], "_oe", sep = "")
  }
}

# okay, now to work that into something for ggplot
TFT.generank.forheatmap.new.forggplot <- data.frame(gene = TFT.generank.forheatmap.new[,1], 
                                                    rank = TFT.generank.forheatmap.new[,2],
                                                    experiment = rep(colnames(TFT.generank.forheatmap.new)[2], 961))
for(i in 1:16){
  tmp <- data.frame(gene = TFT.generank.forheatmap.new[,1], 
                    rank = TFT.generank.forheatmap.new[,(i + 2)],
                    experiment = rep(colnames(TFT.generank.forheatmap.new)[(i + 2)], 961))
  TFT.generank.forheatmap.new.forggplot <- rbind(TFT.generank.forheatmap.new.forggplot, tmp)
}

# now to plot
heatplot <- ggplot(data = TFT.generank.forheatmap.new.forggplot, aes(x = experiment, y = gene, fill = rank)) + 
  geom_tile() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# convert to numeric
TFT.gr.fhm.numeric <- TFT.generank.forheatmap.new[,-1]
for(i in 1:17){
  TFT.gr.fhm.numeric[,i] <- as.numeric(TFT.gr.fhm.numeric[,i])
}

#plot again
heatmap.2(as.matrix(TFT.gr.fhm.numeric), trace = "none", density.info = "none", 
          key = T, col = magma(256))
# try just the zscores directly
heatmap.2(as.matrix(merge.tft.ko.oe[,-1]), trace = "none", density.info = "none",
          key = T, col = viridis(256))


# Work on unique B genes

colnames(uniqueB.new) <- "pr_gene_symbol"
merge.uniqB.ko <- merge(uniqueB.new, TEADko[,c(2,10:18)])
merge.uniqB.ko.oe <- merge(merge.uniqB.ko, TEADoe[,c(2,10:17)], by = "pr_gene_symbol")

# initialize data.frame
uniqB.generank.forheatmap <- data.frame(gene = merge.uniqB.ko.oe[,1], stringsAsFactors = F)
for(i in 1:17){
  print(i)
  if(i <= 9){
    uniqB.generank.forheatmap[,(i+1)] <- order(merge.uniqB.ko.oe[,(i+1)])
    colnames(uniqB.generank.forheatmap)[(i+1)] <- paste(colnames(ko.ordered)[i], "_ko", sep = "")
    
  }else{
    uniqB.generank.forheatmap[,(i+1)] <- order(merge.uniqB.ko.oe[,(i+1)])
    colnames(uniqB.generank.forheatmap)[(i+1)] <- paste(colnames(oe.ordered)[(i-9)], "_oe", sep = "")
  }
}

# make the rank heatmap
uniqB.gr.fhm.numeric <- uniqB.generank.forheatmap[,-1]
rownames(uniqB.gr.fhm.numeric) <- uniqB.generank.forheatmap[,1]
for(i in 1:17){
  uniqB.gr.fhm.numeric[,i] <- as.numeric(uniqB.gr.fhm.numeric[,i])
}

pdf("uniqB_generank_heatmap.pdf")
heatmap.2(as.matrix(uniqB.gr.fhm.numeric), trace = "none", density.info = "none", 
          key = T, col = redblue(100), main = "Unique B cell genes")
dev.off()

# just plot the zscores
rownames(merge.uniqB.ko.oe) <- merge.uniqB.ko.oe[,1]
colnames(merge.uniqB.ko.oe) <- colnames(uniqB.generank.forheatmap)

pdf("uniqB_zscore_heatmap.pdf")
heatmap.2(as.matrix(merge.uniqB.ko.oe[,-1]), trace = "none", density.info = "none", key = T, 
          col = redblue(100), main = "Unique B cell genes, z-scores")
dev.off()

# now unique CD4
colnames(uniqueCD4.new) <- "pr_gene_symbol"

merge.uniqCD4.ko <- merge(uniqueCD4.new, TEADko[,c(2,10:18)])
merge.uniqCD4.ko.oe <- merge(merge.uniqCD4.ko, TEADoe[,c(2,10:17)], by = "pr_gene_symbol")
rownames(merge.uniqCD4.ko.oe) <- merge.uniqCD4.ko.oe[,1]
colnames(merge.uniqCD4.ko.oe) <- colnames(uniqB.generank.forheatmap)


# initialize data.frame
uniqCD4.generank.forheatmap <- data.frame(gene = merge.uniqCD4.ko.oe[,1], stringsAsFactors = F)
for(i in 1:17){
  print(i)
  if(i <= 9){
    uniqCD4.generank.forheatmap[,(i+1)] <- order(merge.uniqCD4.ko.oe[,(i+1)])
    colnames(uniqCD4.generank.forheatmap)[(i+1)] <- paste(colnames(ko.ordered)[i], "_ko", sep = "")
    
  }else{
    uniqCD4.generank.forheatmap[,(i+1)] <- order(merge.uniqCD4.ko.oe[,(i+1)])
    colnames(uniqCD4.generank.forheatmap)[(i+1)] <- paste(colnames(oe.ordered)[(i-9)], "_oe", sep = "")
  }
}

# make the rank heatmap
uniqCD4.gr.fhm.numeric <- uniqCD4.generank.forheatmap[,-1]
rownames(uniqCD4.gr.fhm.numeric) <- uniqCD4.generank.forheatmap[,1]
for(i in 1:17){
  uniqCD4.gr.fhm.numeric[,i] <- as.numeric(uniqCD4.gr.fhm.numeric[,i])
}


pdf("uniqCD4_generank_heatmap.pdf")
heatmap.2(as.matrix(uniqCD4.gr.fhm.numeric), trace = "none", density.info = "none", 
          key = T, col = redblue(100), main = "Unique CD4+ T genes")
dev.off()

heatmap.2(as.matrix(merge.uniqCD4.ko.oe[,-1]), trace = "none", density.info = "none", 
          key = T, col = redblue(100), main = "Unique CD4+ T genes, z-scores")
pdf("uniqCD4_zscore_heatmap.pdf")
heatmap.2(as.matrix(merge.uniqCD4.ko.oe[,-1]), trace = "none", density.info = "none", 
          key = T, col = redblue(100), main = "Unique CD4+ T genes, z-scores", labRow = rep("",154))
dev.off()

# heatmaps of TST and ZMAT5
which(TEADko$pr_gene_symbol == "TST")
#[1] 4017
which(TEADko$pr_gene_symbol == "ZMAT5")
#[1] 10771

TST.heatmapmat <- matrix(nrow = 2, ncol = 9, data = NA)
TST.heatmapmat[1,] <- as.numeric(TEADko[4017,-c(1:9)])
TST.heatmapmat[2,1:3] <- as.numeric(TEADoe[4017,c(10:12)])
TST.heatmapmat[2,5:9] <- as.numeric(TEADoe[4017,c(13:17)])
rownames(TST.heatmapmat) <- c("KD", "OE")
colnames(TST.heatmapmat) <- c("A375", "A549", "HA1E", "HCC515", "HEPG2", "HT29", "MCF7", "PC3", "VCAP")

pdf("TST_zscore_heatmap.pdf")
heatmap.2(TST.heatmapmat, trace = "none", density.info = "none",
          key = T, col = redblue(100), main = "TST z-score", cexRow = 1, na.rm = F)
dev.off()

# and now the same with ZMAT5
ZMAT5.heatmapmat <- matrix(nrow = 2, ncol = 9, data = NA)
ZMAT5.heatmapmat[1,] <- as.numeric(TEADko[10771,-c(1:9)])
ZMAT5.heatmapmat[2,1:3] <- as.numeric(TEADoe[10771,c(10:12)])
ZMAT5.heatmapmat[2,5:9] <- as.numeric(TEADoe[10771,c(13:17)])
rownames(ZMAT5.heatmapmat) <- c("KD", "OE")
colnames(ZMAT5.heatmapmat) <- c("A375", "A549", "HA1E", "HCC515", "HEPG2", "HT29", "MCF7", "PC3", "VCAP")

pdf("ZMAT5_zscore_heatmap.pdf")
heatmap.2(ZMAT5.heatmapmat, trace = "none", density.info = "none",
          key = T, col = redblue(100), main = "ZMAT5 z-score", cexRow = 1, na.rm = F)
dev.off()

# change that up so they are together
TST.ZMAT.mat <- matrix(nrow = 17, ncol = 2, data = NA)
TST.ZMAT.mat[,1] <- as.numeric(c(TST.heatmapmat[1,], TST.heatmapmat[2,1:3], TST.heatmapmat[2,5:9]))
TST.ZMAT.mat[,2] <- as.numeric(c(ZMAT5.heatmapmat[1,], ZMAT5.heatmapmat[2,1:3], ZMAT5.heatmapmat[2,5:9]))
colnames(TST.ZMAT.mat) <- c("TST", "ZMAT5")
rownames(TST.ZMAT.mat) <- c(paste(colnames(ZMAT5.heatmapmat), "KD", sep = " "), 
                            paste(colnames(ZMAT5.heatmapmat)[-4], "OE", sep = " "))

heatmap.2(t(TST.ZMAT.mat), trace = "none", density.info = "none", key = T,
          col = brewer.pal(11, "RdBu"), main = "TST and ZMAT5 z-scores", cexRow = 1)

pdf("TSTandZMAT_zscore_heatmap.pdf")
heatmap.2(t(TST.ZMAT.mat), trace = "none", density.info = "none", key = T,
          col = redblue(100), main = "TST and ZMAT5 z-scores", cexRow = 1, 
          dendrogram = "column", cexCol = 0.85)
dev.off()

heatmap.2(t(TST.ZMAT.mat), trace = "none", density.info = "none", key = T,
          col = plasma(256), main = "TST and ZMAT5 z-scores", cexRow = 1)

# now look at the ranks for those two genes instead of the z-scores
TST.ZMAT.rank <- matrix(nrow = 2, ncol = 17, data = NA)
TST.ZMAT.rank[1,] <- as.numeric(c(ko.ordered[4017,-10], 
                                  oe.ordered[4017,-9]))
TST.ZMAT.rank[2,] <- as.numeric(c(ko.ordered[10771,-10], 
                                  oe.ordered[10771,-9]))
rownames(TST.ZMAT.rank) <- c("TST", "ZMAT5")
colnames(TST.ZMAT.rank) <- rownames(TST.ZMAT.mat)

heatmap.2(TST.ZMAT.rank, trace = "none", density.info = "none", key = T, 
          col = redblue(100), main = "TST and ZMAT5 by rank", cexRow = 1, cexCol = 0.85,
          dendrogram = "column")

pdf("TSTandZMAT_rank_heatmap.pdf")
heatmap.2(TST.ZMAT.rank, trace = "none", density.info = "none", key = T, 
          col = redblue(100), main = "TST and ZMAT5 by rank", cexRow = 1, cexCol = 0.85,
          dendrogram = "column")
dev.off()

png("TSTandZMAT5_rank_heatmap.png", width = 1050, height = 1050, res = 150)
heatmap.2(TST.ZMAT.rank, trace = "none", density.info = "none", key = T, 
          col = redblue(100), main = "TST and ZMAT5 by rank", cexRow = 1, cexCol = 0.85,
          dendrogram = "column")
dev.off()





# okay, now want to move on to the TFT intersect with CD4 and B

uniq.tft.cd4.b <- unique(c(as.character(tft.cd4.new.intersect[,1]), as.character(TFT.B.inter.new[,1])))
uniq.tft.cd4.b <- uniq.tft.cd4.b[which(uniq.tft.cd4.b %in% TEADoe$pr_gene_symbol)]
uniq.tft.cd4.b.heatmapmat <- matrix(nrow = 33, ncol = 17, data = NA)
colnames(uniq.tft.cd4.b.heatmapmat) <- colnames(TST.ZMAT.rank)
rownames(uniq.tft.cd4.b.heatmapmat) <- uniq.tft.cd4.b

for(i in 1:length(uniq.tft.cd4.b)){
  which.row <- which(TEADko$pr_gene_symbol == uniq.tft.cd4.b[i])
  uniq.tft.cd4.b.heatmapmat[i,] <- as.numeric(c(ko.ordered[which.row,-10], 
                                                oe.ordered[which.row,-9]))
}

pdf("Causalgenes_generank.pdf")
heatmap.2(uniq.tft.cd4.b.heatmapmat, trace = "none", density.info = "none", key = T,
          col = redblue(100), main = "Causal genes", cexCol = 0.85)
dev.off()



# okay, now check any of the CD4 and B genes in l1000

uniq.cd4.b <- unique(c(as.character(CD4Tgenes.new[,1]), as.character(Bgenes.new[,1])))
uniq.cd4.b <- uniq.cd4.b[which(uniq.cd4.b %in% TEADoe$pr_gene_symbol)]
uniq.cd4.b.heatmap <- matrix(nrow = 313, ncol = 17, data = NA)
rownames(uniq.cd4.b.heatmap) <- uniq.cd4.b
colnames(uniq.cd4.b.heatmap) <- colnames(uniq.tft.cd4.b.heatmapmat)

for(i in 1:length(uniq.cd4.b)){
  which.row <- which(TEADko$pr_gene_symbol == uniq.cd4.b[i])
  uniq.cd4.b.heatmap[i,] <- as.numeric(c(ko.ordered[which.row,-10], 
                                         oe.ordered[which.row,-9]))
}

pdf("CD4andB_generank.pdf")
heatmap.2(uniq.cd4.b.heatmap, trace = "none", density.info = "none", key = T,
          col = redblue(100), main = "CD4+ T and B cell genes", cexCol = 0.85, labRow =" ")
dev.off()

# then also a z-scores

uniq.cd4.b.heatmap.zscore  <- uniq.cd4.b.heatmap 
for(i in 1:dim(uniq.cd4.b.heatmap.zscore)[1]){
  uniq.cd4.b.heatmap.zscore[i,] <-  c(as.numeric(
                                        TEADko[which(TEADko$pr_gene_symbol == rownames(uniq.cd4.b.heatmap.zscore)[i]),-c(1:9)]),
                                      as.numeric(
                                        TEADoe[which(TEADoe$pr_gene_symbol == rownames(uniq.cd4.b.heatmap.zscore)[i]),-c(1:9)]))
}

pdf("CD4andB_zscore.pdf")
heatmap.2(uniq.cd4.b.heatmap.zscore, trace = "none", density.info = "none", key = T,
          col = redblue(100), main = "CD4+ T and B cell genes", cexCol = 0.85, labRow =" ", Colv = F, dendrogram = "row")
dev.off()



# want to figure out what genes had consistent patterns in the cd4 and B genes 
# TST is what we want, higher zscores in OE and lower KD

list.tstlike <- c()
for(i in 1:313){
  if(sum(uniq.cd4.b.heatmap.zscore[i,10:17] >= 0) >= 6){
    if(sum(uniq.cd4.b.heatmap.zscore[i, 1:9] < 0) >= 7){
      list.tstlike <- c(list.tstlike, rownames(uniq.cd4.b.heatmap.zscore)[i])
    }
  }
}
# if i say that all 8 oes must be the same and all 9 kos must be the same, don't even get tst 
# (one oe is slightly less than 0), if I fudge it so that one ko and oe can be off, only get tst, and if
# I make it so that two of each can be off I get 7 results


pdf("CD4andB_samedirections_zscore.pdf")
heatmap.2(uniq.cd4.b.heatmap.zscore[c(list.tstlike, oppo.list),], trace = "none", density.info = "none", key = T,
          col = redblue(100), main = "CD4+ T and B cell genes", cexCol = 0.85, Colv = F, dendrogram = "row")
dev.off()

# can also check if there are any that are the other way
oppo.list <- c()
for(i in 1:313){
  if(sum(uniq.cd4.b.heatmap.zscore[i,10:17] < 0) >= 6){
    if(sum(uniq.cd4.b.heatmap.zscore[i, 1:9] > 0) >= 7){
      oppo.list <- c(oppo.list, rownames(uniq.cd4.b.heatmap.zscore)[i])
    }
  }
}

# got 4 hits there
pdf("CD4andB_sameandoppodirection_zscore.pdf")
heatmap.2(uniq.cd4.b.heatmap.zscore[c(list.tstlike,oppo.list)], trace = "none", density.info = "none", key = T,
          col = redblue(100), main = "CD4+ T and B cell genes", cexCol = 0.85, Colv = F, dendrogram = "row")
dev.off()

# lets make a matrix that says if the z scores are up or down for each gene, and then have
# counts for each type of experiement which are up or down

cd4.updown <- as.data.frame(matrix(NA, ncol = 21, nrow = length(CD4Tgenes.new[(which(CD4Tgenes.new[,1] %in% TEADko$pr_gene_symbol)),1])), row.names = CD4Tgenes.new[(which(CD4Tgenes.new[,1] %in% TEADko$pr_gene_symbol)),1])
colnames(cd4.updown) <- c(colnames(uniq.cd4.b.heatmap.zscore), "KO_up", "KO_down", "OE_up", "OE_down")
for(i in 1:266){
  ko.row <- TEADko[which(TEADko$pr_gene_symbol == rownames(cd4.updown)[i]), -c(1:9)] 
  oe.row <- TEADoe[which(TEADoe$pr_gene_symbol == rownames(cd4.updown)[i]), -c(1:9)] 
  for(j in 1:9){
    if(ko.row[j] < 0){
      cd4.updown[i,j] <- -1
    }else if(ko.row[j] > 0){
      cd4.updown[i,j] <- 1
    }else{cd4.updown[i,j] <- 0}
  }
  for(k in 1:8){
    if(oe.row[k] < 0){
      cd4.updown[i, (k + 9)] <- -1
    }else if(oe.row[k] > 0){
      cd4.updown[i, (k + 9)] <- 1
    }else{cd4.updown[i, (k + 9)] <- 0}
  }
}

cd4.updown$KO_up <- apply(cd4.updown[,1:9], 1, function(x) sum(x == 1))
cd4.updown$KO_down <- apply(cd4.updown[,1:9], 1, function(x) sum(x == -1))
cd4.updown$OE_up <- apply(cd4.updown[,10:17], 1, function(x) sum(x == 1))
cd4.updown$OE_down <- apply(cd4.updown[,10:17], 1, function(x) sum(x == -1))

# which are all up in KO and down in OE
for(i in 1:266){
  if(sum(cd4.updown$KO_up[i], cd4.updown$OE_down[i]) == 17){
    print(rownames(cd4.updown)[i])
  }
}
#none


#which are all down in KO and up in OE
for(i in 1:266){
  if(sum(cd4.updown$KO_down[i], cd4.updown$OE_up[i]) == 17){
    print(rownames(cd4.updown)[i])
  }
}
#none

# off by one?
for(i in 1:266){
if(sum(cd4.updown$KO_up[i], cd4.updown$OE_down[i]) >= 16){
  print(c(rownames(cd4.updown)[i], "  kUoD"))
}else if(sum(cd4.updown$KO_down[i], cd4.updown$OE_up[i]) >= 16){
  print(c(rownames(cd4.updown)[i], "  kDoU"))
}}
# just TST

# off by two?
for(i in 1:266){
  if(sum(cd4.updown$KO_up[i], cd4.updown$OE_down[i]) >= 15){
    print(c(rownames(cd4.updown)[i], "kUoD"))
  }else if(sum(cd4.updown$KO_down[i], cd4.updown$OE_up[i]) >= 15){
    print(c(rownames(cd4.updown)[i], "kDoU"))
  }}
# [1] "TMEM109" "kUoD" 
# [1] "TM9SF2" "kUoD"
# [1] "TST"    "kDoU"

## do the same for the B list 
b.updown <- as.data.frame(matrix(NA, ncol = 21, nrow = length(Bgenes.new[(which(Bgenes.new[,1] %in% TEADko$pr_gene_symbol)),1])), row.names = Bgenes.new[(which(Bgenes.new[,1] %in% TEADko$pr_gene_symbol)),1])
colnames(b.updown) <- c(colnames(uniq.cd4.b.heatmap.zscore), "KO_up", "KO_down", "OE_up", "OE_down")
for(i in 1:174){
  ko.row <- TEADko[which(TEADko$pr_gene_symbol == rownames(b.updown)[i]), -c(1:9)] 
  oe.row <- TEADoe[which(TEADoe$pr_gene_symbol == rownames(b.updown)[i]), -c(1:9)] 
  for(j in 1:9){
    if(ko.row[j] < 0){
      b.updown[i,j] <- -1
    }else if(ko.row[j] > 0){
      b.updown[i,j] <- 1
    }else{b.updown[i,j] <- 0}
  }
  for(k in 1:8){
    if(oe.row[k] < 0){
      b.updown[i, (k + 9)] <- -1
    }else if(oe.row[k] > 0){
      b.updown[i, (k + 9)] <- 1
    }else{b.updown[i, (k + 9)] <- 0}
  }
}

b.updown$KO_up <- apply(b.updown[,1:9], 1, function(x) sum(x == 1))
b.updown$KO_down <- apply(b.updown[,1:9], 1, function(x) sum(x == -1))
b.updown$OE_up <- apply(b.updown[,10:17], 1, function(x) sum(x == 1))
b.updown$OE_down <- apply(b.updown[,10:17], 1, function(x) sum(x == -1))

# which are all up in one and down in the other
for(i in 1:174){
  if(sum(b.updown$KO_up[i], b.updown$OE_down[i]) == 17){
    print(c(rownames(b.updown)[i], "  kUoD"))
  }else if(sum(b.updown$KO_down[i], b.updown$OE_up[i]) == 17){
    print(c(rownames(b.updown)[i], "  kDoU"))
  }}
#none

# off by one?
for(i in 1:174){
  if(sum(b.updown$KO_up[i], b.updown$OE_down[i]) >= 16){
    print(c(rownames(b.updown)[i], "  kUoD"))
  }else if(sum(b.updown$KO_down[i], b.updown$OE_up[i]) >= 16){
    print(c(rownames(b.updown)[i], "  kDoU"))
  }}
# none

# off by two?
for(i in 1:174){
  if(sum(b.updown$KO_up[i], b.updown$OE_down[i]) >= 15){
    print(c(rownames(b.updown)[i], "kUoD"))
  }else if(sum(b.updown$KO_down[i], b.updown$OE_up[i]) >= 15){
    print(c(rownames(b.updown)[i], "kDoU"))
  }}
#none

for(i in 1:174){
  if(sum(b.updown$KO_up[i], b.updown$OE_down[i]) >= 14){
    print(c(rownames(b.updown)[i], "kUoD"))
  }else if(sum(b.updown$KO_down[i], b.updown$OE_up[i]) >= 14){
    print(c(rownames(b.updown)[i], "kDoU"))
  }}
# finally, 2 
# [1] "RPL5" "kDoU"
# [1] "SATB1" "kDoU" 


# combine the two 
cd4.b.updown <- rbind(cd4.updown, b.updown[-which(rownames(b.updown) %in% rownames(cd4.updown)),])

for(i in 1:313){
  if(sum(cd4.b.updown$KO_up[i], cd4.b.updown$OE_down[i]) >= 14){
    print(c(rownames(cd4.b.updown)[i], "kUoD"))
  }else if(sum(cd4.b.updown$KO_down[i], cd4.b.updown$OE_up[i]) >= 14){
    print(c(rownames(cd4.b.updown)[i], "kDoU"))
  }
}

# save this table 
write.table(cd4.b.updown, file = "CD4andB.zscoredirections.txt", quote = F)
write.table(uniq.cd4.b.heatmap.zscore, file = "CD4andB.actualzscores.txt", quote = F)


# create 2 files: one with zscores, one with ranks, for each gene in each cell line, 
# tell whether in TFT list, CD4 list, B list, keep cell lines together

# combine all the tft, cd4, and be genes
#all.tft.cd4.b <- unique(unique(TEADtft[,1], CD4Tgenes.new[,1]), Bgenes.new[,1])
# thats just the full TFT list

zscores.TFT.CD4.B <- as.data.frame(matrix(NA, ncol = 20, nrow = 313))
rownames(zscores.TFT.CD4.B) <- rownames(cd4.b.updown)
colnames(zscores.TFT.CD4.B) <- c("A375_KD", "A375_OE", "A549_KD", "A549_OE", "HA1E_KD", "HA1E_OE", 
                                "HEPG2_KD", "HEPG2_OE", "HT29_KD", "HT29_OE", "MCF7_KD", "MCF7_OE", 
                                "PC3_KD", "PC3_OE", "VCAP_KD", "VCAP_OE", "HCC515_KD", "is_TFT", "CD4", "B")
for(i in 1:313){
  check <- rownames(zscores.TFT.CD4.B)[i]
  if(check %in% TEADtft[,1]){
    zscores.TFT.CD4.B$is_TFT[i] <- 1
  }else{zscores.TFT.CD4.B$is_TFT[i] <- 0}
  
  if(check %in% CD4Tgenes.new[,1]){
    zscores.TFT.CD4.B$CD4[i] <- 1
  }else{zscores.TFT.CD4.B$CD4[i] <- 0}
  
  if(check %in% Bgenes.new[,1]){
    zscores.TFT.CD4.B$B[i] <- 1
  }else{zscores.TFT.CD4.B$B[i] <- 0}
}

# can copy that now for the ranks as well
ranks.TFT.CD4.B <- zscores.TFT.CD4.B

# and now I can fill those two in 
zscores.TFT.CD4.B[1:17] <- uniq.cd4.b.heatmap.zscore[,c(1,10,2,11,3,12,5,13,6,14,7,15,8,16,9,17,4)]
ranks.TFT.CD4.B[1:17] <- uniq.cd4.b.heatmap[,c(1,10,2,11,3,12,5,13,6,14,7,15,8,16,9,17,4)]

zscores.TFT.CD4.B <- data.frame(Gene = rownames(zscores.TFT.CD4.B), zscores.TFT.CD4.B)
ranks.TFT.CD4.B <- data.frame(Gene = rownames(ranks.TFT.CD4.B), ranks.TFT.CD4.B)

# save those files
write.table(zscores.TFT.CD4.B, "zscores.CD4.B.txt", quote = F, row.names = F)
write.table(ranks.TFT.CD4.B, "ranks.CD4.B.txt", quote = F, row.names = F)

# save workspace
save(list = ls(all.names = T), file = "TEAD2_prioritized.RData")
