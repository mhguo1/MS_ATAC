###### Pathway and PPI analyses
## last update: 2/12/2021

### libraries
library(ggplot2)
library(ggrepel)
library(VennDiagram)


### working directory
setwd("/Users/np911/Dropbox (Partners HealthCare)/Projects/MS-chromatin_fine_mapping/results/Fine-mapping")

### read data | canonical pathways results

## CD4 
cd4.data <- read.table("gsea/CD4_T/cp.enrichment", header=TRUE, stringsAsFactors = FALSE, sep = "\t")
dim(cd4.data)
names(cd4.data)[c(2,3,4,6)] <- c("CD4_pval", "CD4_qval", "CD4_numSharedGenes", "CD4_sharedGenes")
# dataframe for join
cd4.sub <- cd4.data[ , c(1, 3, 4)]
cd4.sub$category <- "CD4 T"
names(cd4.sub)[c(2,3)] <- c("qval", "numSharedGenes")

# B cells
b.data <- read.table("gsea/B/cp.enrichment", header=TRUE, stringsAsFactors = FALSE, sep = "\t")
dim(b.data)
names(b.data)[c(2,3,4,6)] <- c("B_pval", "B_qval", "B_numSharedGenes", "B_sharedGenes")
# dataframe for join
b.sub <- b.data[ , c(1, 3, 4)]
b.sub$category <- "B"
names(b.sub)[c(2,3)] <- c("qval", "numSharedGenes")

# CD4 B common
cd4_b.data <- read.table("gsea/CD4_B_common/cp.enrichment", header=TRUE, stringsAsFactors = FALSE, sep = "\t")
dim(cd4_b.data)
names(cd4_b.data)[c(2,3,4,6)] <- c("cd4_b_pval", "cd4_b_qval", "cd4_b_numSharedGenes", "cd4_b_sharedGenes")
# dataframe for join
cd4_b.sub <- cd4_b.data[ , c(1, 3, 4)]
cd4_b.sub$category <- "Common"
names(cd4_b.sub)[c(2,3)] <- c("qval", "numSharedGenes")


# CD4 unique
cd4Un.data <- read.table("gsea/CD4_unique/cp.enrichment", header=TRUE, stringsAsFactors = FALSE, sep = "\t")
dim(cd4Un.data)
names(cd4Un.data)[c(2,3,4,6)] <- c("cd4Un_pval", "cd4Un_qval", "cd4Un_numSharedGenes", "cd4Un_sharedGenes")
# dataframe for join
cd4Un.sub <- cd4Un.data[ , c(1, 3, 4)]
cd4Un.sub$category <- "CD4 T Unique"
names(cd4Un.sub)[c(2,3)] <- c("qval", "numSharedGenes")


# B unique
bUn.data <- read.table("gsea/B_unique/cp.enrichment", header=TRUE, stringsAsFactors = FALSE, sep = "\t")
dim(bUn.data)
names(bUn.data)[c(2,3,4,6)] <- c("bUn_pval", "bUn_qval", "bUn_numSharedGenes", "bUn_sharedGenes")
# dataframe for join
bUn.sub <- bUn.data[ , c(1, 3, 4)]
bUn.sub$category <- "B Unique"
names(bUn.sub)[c(2,3)] <- c("qval", "numSharedGenes")


### create a joined dataframe
all.sub <-rbind(cd4.sub, b.sub, cd4_b.sub, cd4Un.sub, bUn.sub)
## create a column with factorized -log10(P)
all.sub$FDR <- -log10(all.sub$qval)
max(all.sub$FDR)
all.sub$FDR<-cut(all.sub$FDR,breaks=c(0, 1.30103, 2, 3, 4, 5, max(all.sub$FDR)),include.lowest=TRUE,label=c(" > 5e-2","5e-2 to 1e-2", "1e-2 to 1e-3","1e-3 to 1e-4","1e-4 to 1e-5"," < 1e-6")) 
# reorder the "levels" of category
all.sub$category <- factor(all.sub$category, levels = c("CD4 T", "B", "Common", "CD4 T Unique", "B Unique"))
## keep only the genesets with at least one qvalue<0.05
all.sub <- all.sub[ all.sub$geneSetName %in% unique(all.sub$geneSetName[all.sub$qval<0.05]), ] 



### heatmap
#create heatmap using rescaled values
ggplot(all.sub, aes(category,geneSetName)) +
  geom_tile(aes(fill = FDR)) +
  scale_fill_manual(values = c("white", "grey50", "grey40", "grey30", "grey20", "grey10"),name="FDR") +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 90), axis.ticks.y = element_blank()) +
  labs(x="Prioritized genes", y="Canonical Pathways")

### Venn diagram of CD4 T and B cells
#cd4_genes <- read.table("CD4_T_cells_genes.txt", header=FALSE, stringsAsFactors = FALSE)
#b_genes <- read.table("B_cells_genes.txt", header=FALSE, stringsAsFactors = FALSE)
cd4_genes <- read.table("CD4_prioritized_genes.txt", header=FALSE, stringsAsFactors = FALSE)
b_genes <- read.table("B_prioritized_genes.txt", header=FALSE, stringsAsFactors = FALSE)

# create a list 
venn_list <- list(CD4=cd4_genes$V1, B=b_genes$V1)
# venn diagram
grid.newpage()
draw.pairwise.venn(area1 = nrow(cd4_genes),
                                area2      = nrow(b_genes),
                                cross.area = sum(cd4_genes$V1 %in% b_genes$V1),
                                fill            = c('#E69F00', '#999999'),
                   cex = 0
                                )

### Pathway comparisons
## CD4_unique vs. B_unique
# merge
joined <- merge(cd4Un.data, bUn.data)
dim(joined)
# create a category for the plots. Note that there is no pathway with FDR < 0.05 in either cells.
joined$category <- ifelse(joined$cd4Un_qval<0.05, "CD4 T", "None" )
joined$category[joined$bUn_qval<0.05] <- "B"
joined$category <- as.factor(joined$category)

#find the minimum qvalue 
min_all <- min( min(joined$cd4Un_qval), min(joined$bUn_qval))
# plot with labels for top pathways
ggplot(joined, aes(x=-log10(cd4Un_qval), y=-log10(bUn_qval),  label=geneSetName, fill=category)) + 
  geom_point(shape=21, aes(size=numAllGenesInSet)) +
  labs(size="Genes in pathway") +
  scale_fill_manual(values=c('#999999','#E69F00', 'lightgrey')) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", 
             color = "red", size=0.2) +
  geom_vline(xintercept=-log10(0.05), linetype="dashed", 
             color = "red", size=0.2) +
  geom_abline(intercept = 0, slope = 1, size=0.2) +
  xlim(0, (-log10(min_all)+ 1)) +
  ylim(0, (-log10(min_all)+ 1) ) +
  theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill = "white", colour = "black",
                                                                                size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "grey"), ) +
  xlab("Genes unique in CD4  cells") +
  ylab("Genes unique in B cells") +
  ggtitle("Canonical pathway enrichment") +
  geom_label_repel(data         = subset(joined, -log10(bUn_qval) > 3.5 | -log10(cd4Un_qval) > 4 ),
                   size          = 2,
                   box.padding   = 0.3,
                   point.padding = 0.8,
                   force         = 10,
                   segment.size  = 0.5,
                   segment.color = "grey50",
                   direction     = "y",
                   show.legend=F) 

# plot without labels
ggplot(joined, aes(x=-log10(cd4Un_qval), y=-log10(bUn_qval),  label=geneSetName, fill=category)) + 
  geom_point(shape=21, aes(size=numAllGenesInSet)) +
  labs(size="Genes in pathway") +
  scale_fill_manual(values=c('#999999','#E69F00', 'lightgrey')) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", 
             color = "red", size=0.2) +
  geom_vline(xintercept=-log10(0.05), linetype="dashed", 
             color = "red", size=0.2) +
  geom_abline(intercept = 0, slope = 1, size=0.2) +
  xlim(0, (-log10(min_all)+ 1)) +
  ylim(0, (-log10(min_all)+ 1) ) +
  theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill = "white", colour = "black",
                                                                                size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "grey"), ) +
  xlab("Genes unique in CD4  cells") +
  ylab("Genes unique in B cells") +
  ggtitle("Canonical pathway enrichment")


  
## create joint dataframe for Sup Table
tmp1 <- merge(cd4.data, b.data)
nrow(tmp1)
tmp2 <- merge(tmp1, cd4_b.data)
nrow(tmp2)
tmp3 <- merge(tmp2, cd4Un.data)
nrow(tmp3)
all <- merge(tmp3, bUn.data)
nrow(all)
# save the table
write.table(all, "gsea/all_pathways.txt", row.names = FALSE, quote = FALSE, sep = ",")
