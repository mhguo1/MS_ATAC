---
title: "TEAD2_figures updated,prioritized list"
author: "Nikolaos Patsopoulos; edited Brenna LaBarre"
date: "02/18/2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
# working directory
require("knitr")
opts_knit$set(root.dir = "~/Dropbox (Partners HealthCare)/MS-chromatin_fine_mapping/")
```

## Load libraries  


```{r Environment setup}
## Libraries
library(ggplot2)
library(ggrepel)



```

## TFT plots
Organize the data

```{r TFT, echo=FALSE}
## CD4
cd4.tft.gsea <- read.table("results/Fine-mapping/gsea/CD4_T/gtrd.enrichment", sep="\t", header=TRUE, stringsAsFactors = FALSE)
# keep only few columns
cd4.tft.gsea <- cd4.tft.gsea[ , c(1, 2, 3, 4, 5)]
# add percent shared genes
cd4.tft.gsea$perSharedGenes <- cd4.tft.gsea$numSharedGenes / cd4.tft.gsea$numAllGenesInSet
# add value with category
cd4.tft.gsea$category <- "CD4+ T"
# add variable with color 
cd4.tft.gsea$cols <- c("lightgrey", "red")[(cd4.tft.gsea$geneSetName=="TEAD2_TARGET_GENES") + 1]  
# subset only the qvalue<0.05
cd4.tft.gsea.sub <- cd4.tft.gsea[ cd4.tft.gsea$qval<0.05, ]
#barplot(-log10(cd4.tft[cd4.tft$qval<0.05, 3]), col=cols)

## B
b.tft.gsea <- read.table("results/Fine-mapping/gsea/B/gtrd.enrichment", sep="\t", header=TRUE, stringsAsFactors = FALSE)
# keep only few columns
b.tft.gsea <- b.tft.gsea[ , c(1, 2, 3, 4, 5)]
# add percent shared genes
b.tft.gsea$perSharedGenes <- b.tft.gsea$numSharedGenes / b.tft.gsea$numAllGenesInSet
# add value with category
b.tft.gsea$category <- "B"
# add variable with color 
b.tft.gsea$cols <- c("lightgrey", "red")[(b.tft.gsea$geneSetName=="TEAD2_TARGET_GENES") + 1]  
# subset only the qvalue<0.05
b.tft.gsea.sub <- b.tft.gsea[ b.tft.gsea$qval<0.05, ]

## Common
common.tft.gsea <- read.table("results/Fine-mapping/gsea/CD4_B_common/gtrd.enrichment", sep="\t", header=TRUE, stringsAsFactors = FALSE)
# keep only few columns
common.tft.gsea <- common.tft.gsea[ , c(1, 2, 3, 4, 5)]
# add percent shared genes
common.tft.gsea$perSharedGenes <- common.tft.gsea$numSharedGenes / common.tft.gsea$numAllGenesInSet
# add value with category
common.tft.gsea$category <- "Common"
# add variable with color 
common.tft.gsea$cols <- c("lightgrey", "red")[(common.tft.gsea$geneSetName=="TEAD2_TARGET_GENES") + 1]  
# subset only the qvalue<0.05
common.tft.gsea.sub <- common.tft.gsea[ common.tft.gsea$qval<0.05, ]

## CD4 unique
cd4Un.tft.gsea <- read.table("results/Fine-mapping/gsea/CD4_unique/gtrd.enrichment", sep="\t", header=TRUE, stringsAsFactors = FALSE)
# keep only few columns
cd4Un.tft.gsea <- cd4Un.tft.gsea[ , c(1, 2, 3, 4, 5)]
# add percent shared genes
cd4Un.tft.gsea$perSharedGenes <- cd4Un.tft.gsea$numSharedGenes / cd4Un.tft.gsea$numAllGenesInSet
# add value with category
cd4Un.tft.gsea$category <- "CD4+ T Unique"
# add variable with color 
cd4Un.tft.gsea$cols <- c("lightgrey", "red")[(cd4Un.tft.gsea$geneSetName=="TEAD2_TARGET_GENES") + 1]  
# subset only the qvalue<0.05
cd4Un.tft.gsea.sub <- cd4Un.tft.gsea[ cd4Un.tft.gsea$qval<0.05, ]

## B unique
bUn.tft.gsea <- read.table("results/Fine-mapping/gsea/B_unique/gtrd.enrichment", sep="\t", header=TRUE, stringsAsFactors = FALSE)
# keep only few columns
bUn.tft.gsea <- bUn.tft.gsea[ , c(1, 2, 3, 4, 5)]
# add percent shared genes
bUn.tft.gsea$perSharedGenes <- bUn.tft.gsea$numSharedGenes / bUn.tft.gsea$numAllGenesInSet
# add value with category
bUn.tft.gsea$category <- "B Unique"
# add variable with color 
bUn.tft.gsea$cols <- c("lightgrey", "red")[(bUn.tft.gsea$geneSetName=="TEAD2_TARGET_GENES") + 1]  
# subset only the qvalue<0.05
bUn.tft.gsea.sub <- bUn.tft.gsea[ bUn.tft.gsea$qval<0.05, ]


## join them
all.tft.gsea <- rbind(cd4.tft.gsea, b.tft.gsea, common.tft.gsea, cd4Un.tft.gsea, bUn.tft.gsea)
## order based on qvalue
all.tft.gsea <- all.tft.gsea[order(all.tft.gsea$qval), ]

```

Plots
Note: facet doesn't sort properly the bars. We will use grid.arrange instead.
```{r TFT plot, echo=FALSE}
element_textbox <- function(...) {
  el <- element_text(...)
  class(el) <- c("element_textbox", class(el))
  el
}

element_grob.element_textbox <- function(element, ...) {
  text_grob <- NextMethod()
  rect_grob <- element_grob(calc_element("strip.background", theme_bw()))
  
  ggplot2:::absoluteGrob(
    grid::gList(
      element_grob(calc_element("strip.background", theme_bw())),
      text_grob
    ),
    height = grid::grobHeight(text_grob), 
    width = grid::unit(1, "npc")
  )
}


#all.tft$cols <- c("lightgrey", "red")[(all.tft$geneSetName=="TEAD2_TARGET_GENES") + 1]  
#ggplot(all.tft, aes(reorder(geneSetName, -log10(qval)), -log10(qval))) +
#  geom_col(col=all.tft$cols) +
#  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#  facet_wrap(~category, nrow = 5)

# cd4
p1 <- ggplot(cd4.tft.gsea.sub, aes(reorder(geneSetName, qval), -log10(qval))) +
  geom_col(col=cd4.tft.gsea.sub$cols) +
 theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)) ) +
  labs(y="-log10(FDR)", x="", title="CD4+ T genes") 
# b 
p2 <- ggplot(b.tft.gsea.sub, aes(reorder(geneSetName, qval), -log10(qval))) +
  geom_col(col=b.tft.gsea.sub$cols) +
 theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)) ) +
  labs(y="-log10(FDR)", x="", title="B genes") 
# Common
p3 <- ggplot(common.tft.gsea.sub, aes(reorder(geneSetName, qval), -log10(qval))) +
  geom_col(col=common.tft.gsea.sub$cols) +
 theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)) ) +
  labs(y="-log10(FDR)", x="", title="Common genes") 


# CD4+ T unique
p4 <- ggplot(cd4Un.tft.gsea.sub, aes(reorder(geneSetName, qval), -log10(qval))) +
  geom_col(col=cd4Un.tft.gsea.sub$cols) +
 theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)) ) +
  labs(y="-log10(FDR)", x="", title="Unique CD4+ T genes") 

# B unique
p5 <- ggplot(bUn.tft.gsea.sub, aes(reorder(geneSetName, qval), -log10(qval))) +
  geom_col(col=bUn.tft.gsea.sub$cols) +
 theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)) ) +
  labs(y="-log10(FDR)", x="", title="Unique B genes") 


```

```{r saveimages}
ggsave("cd4.enrichmentbarplot.pdf", p1, device = "pdf", path = "~/Desktop")
ggsave("B.enrichmentbarplot.pdf", p2, device = "pdf", path = "~/Desktop")
ggsave("common.enrichmentbarplot.pdf", p3, device = "pdf", path = "~/Desktop")
ggsave("cd4Un.enrichmentbarplot.pdf", p4, device = "pdf", path = "~/Desktop")
ggsave("BUn.enrichmentbarplot.pdf", p5, device = "pdf", path = "~/Desktop")

```


```{r otherplot}
# change up to a different style of plots 
# x-axis is rank of enrichment and the y is the log10 of pvalue ( same a above, just not bars)
# make a vector with the labels
#p6.labels <- reorder(cd4.tft.gsea.sub$geneSetName, cd4.tft.gsea.sub$qval)
# remove  the "target_genes"
p6.labels <- gsub("_TARGET_GENES", "", cd4.tft.gsea.sub$geneSetName)
# only keep the TEAD2 lable
which(p6.labels == "TEAD2")
# it's 16
p6.labels[c(1:15, 17:length(p6.labels))] <- ""

cd4.tft.gsea.sub$label <- p6.labels

p6 <- ggplot(cd4.tft.gsea.sub, aes(reorder(geneSetName, qval), -log10(qval), label = label)) +
  geom_point(aes(color = reorder(geneSetName, qval))) + geom_label_repel() +
 theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)) ) +
  labs(y="-log10(FDR)", x="", title="CD4+ T genes") + scale_color_viridis(discrete = T)


# B

b.tft.gsea.sub$label <- gsub("_TARGET_GENES", "", b.tft.gsea.sub$geneSetName)
which(b.tft.gsea.sub$label == "TEAD2") #34
b.tft.gsea.sub$label[c(1:33, 35:length(b.tft.gsea.sub$label))] <- ""

p7 <- ggplot(b.tft.gsea.sub, aes(reorder(geneSetName, qval), -log10(qval), label = label)) +
  geom_point(aes(color = reorder(geneSetName, qval))) + geom_label_repel() +
 theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)) ) +
  labs(y="-log10(FDR)", x="", title="B genes") + scale_color_viridis(discrete = T)

# common
common.tft.gsea.sub$label <- gsub("_TARGET_GENES", "", common.tft.gsea.sub$geneSetName)
which(common.tft.gsea.sub$label == "TEAD2") #8 
common.tft.gsea.sub$label[c(1:7, 9:length(common.tft.gsea.sub$label))] <- ""

p8 <- ggplot(common.tft.gsea.sub, aes(reorder(geneSetName, qval), -log10(qval), label = label)) +
  geom_point(aes(color = reorder(geneSetName, qval))) + geom_label_repel() +
 theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)) ) +
  labs(y="-log10(FDR)", x="", title="Common genes") + scale_color_viridis(discrete = T)


# cd4U
cd4Un.tft.gsea.sub$label <- gsub("_TARGET_GENES", "", cd4Un.tft.gsea.sub$geneSetName)
which(cd4Un.tft.gsea.sub$label == "TEAD2") #70
cd4Un.tft.gsea.sub$label[c(1:69, 71:length(cd4Un.tft.gsea.sub$label))] <- ""

p9 <- ggplot(cd4Un.tft.gsea.sub, aes(reorder(geneSetName, qval), -log10(qval), label = label)) +
  geom_point(aes(color = reorder(geneSetName, qval))) + geom_label_repel() +
 theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)) ) +
  labs(y="-log10(FDR)", x="", title="Unique CD4+ T genes") + scale_color_viridis(discrete = T)


# BU
bUn.tft.gsea.sub$label <- gsub("_TARGET_GENES", "", bUn.tft.gsea.sub$geneSetName)
which(bUn.tft.gsea.sub$label == "TEAD2") # not there, only 30 total genes
bUn.tft.gsea.sub$label <- ""

p10 <- ggplot(bUn.tft.gsea.sub, aes(reorder(geneSetName, qval), -log10(qval), label = label)) +
  geom_point(aes(color = reorder(geneSetName, qval))) + geom_label_repel() +
 theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)) ) +
  labs(y="-log10(FDR)", x="") + scale_color_viridis(discrete = T)

# then save those
ggsave("cd4.enrichmentlabeleddots.pdf", p6, device = "pdf", path = "~/Desktop")
ggsave("B.enrichmentlabeleddots.pdf", p7, device = "pdf", path = "~/Desktop")
ggsave("common.enrichmentlabeleddots.pdf", p8, device = "pdf", path = "~/Desktop")
ggsave("cd4Un.enrichmentlabeleddots.pdf", p9, device = "pdf", path = "~/Desktop")
ggsave("BUn.enrichmentlabeleddots.pdf", p10, device = "pdf", path = "~/Desktop")


## Remake them using the whole list, not the subset; remove title; 

#cd4
# remove  the "target_genes"
p11.labels <- gsub("_TARGET_GENES", "", cd4.tft.gsea$geneSetName)
# only keep the TEAD2 lable
which(p11.labels == "TEAD2")
# it's 16
p11.labels[c(1:15, 17:length(p11.labels))] <- ""

cd4.tft.gsea$label <- p11.labels

p11 <- ggplot(cd4.tft.gsea, aes(reorder(geneSetName, qval), -log10(qval), label = label)) +
  geom_point(color = "black", fill = "white", shape = 21) + geom_label_repel(min.segment.length = 0, point.padding = 1) + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray"), panel.grid.major.x = element_blank(), axis.line.y = element_line(color = "black")) + labs(y="-log10(FDR)", x="") + expand_limits(x = -5) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")

  
#B
b.tft.gsea$label <- gsub("_TARGET_GENES", "", b.tft.gsea$geneSetName)
which(b.tft.gsea$label == "TEAD2") #34
b.tft.gsea$label[c(1:33, 35:length(b.tft.gsea$label))] <- ""
  
p12 <- ggplot(b.tft.gsea, aes(reorder(geneSetName, qval), -log10(qval), label = label)) +
  geom_point(color = "black", fill = "white", shape = 21) + geom_label_repel(min.segment.length = 0, point.padding = 1, nudge_y = 0.5) + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray"), panel.grid.major.x = element_blank(), axis.line.y = element_line(color = "black")) + labs(y="-log10(FDR)", x="") + expand_limits(x = -5) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")

#common
common.tft.gsea$label <- gsub("_TARGET_GENES", "", common.tft.gsea$geneSetName)
which(common.tft.gsea$label == "TEAD2") #8
common.tft.gsea$label[c(1:7, 9:length(common.tft.gsea$label))] <- ""

p13 <- ggplot(common.tft.gsea, aes(reorder(geneSetName, qval), -log10(qval), label = label)) +
  geom_point(color = "black", fill = "white", shape = 21) + geom_label_repel(min.segment.length = 0, point.padding = 1) + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray"), panel.grid.major.x = element_blank(), axis.line.y = element_line(color = "black")) + labs(y="-log10(FDR)", x="") + expand_limits(x = -8) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")

#cd4un
cd4Un.tft.gsea$label <- gsub("_TARGET_GENES", "", cd4Un.tft.gsea$geneSetName)
which(cd4Un.tft.gsea$label == "TEAD2") #70
cd4Un.tft.gsea$label[c(1:69, 71:length(cd4Un.tft.gsea$label))] <- ""

p14 <- ggplot(cd4Un.tft.gsea, aes(reorder(geneSetName, qval), -log10(qval), label = label)) +
  geom_point(color = "black", fill = "white", shape = 21) + geom_label_repel(min.segment.length = 0, point.padding = 1) + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray"), panel.grid.major.x = element_blank(), axis.line.y = element_line(color = "black")) + labs(y="-log10(FDR)", x="") + expand_limits(x = -5) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")

  
#bun

bUn.tft.gsea$label <- gsub("_TARGET_GENES", "", bUn.tft.gsea$geneSetName)
which(bUn.tft.gsea$label == "TEAD2") #316
bUn.tft.gsea$label[c(1:315, 317:length(bUn.tft.gsea$label))] <- ""

p15 <- ggplot(bUn.tft.gsea, aes(reorder(geneSetName, qval), -log10(qval), label = label)) +
  geom_point(color = "black", fill = "white", shape = 21) + geom_label_repel(min.segment.length = 0, point.padding = 1, nudge_y = 0.5) + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray"), panel.grid.major.x = element_blank(), axis.line.y = element_line(color = "black")) + labs(y="-log10(FDR)", x="") + expand_limits(x = -8) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")

# And now save those
ggsave("cd4.enrichmentlabeleddots.FullTFT.pdf", p11, device = "pdf", path = "/Users/brenna/Dropbox (Partners HealthCare)/MS-chromatin_fine_mapping/results/Fine-mapping/TEAD2_prioritized_Brenna/")
ggsave("B.enrichmentlabeleddots.FullTFT.pdf", p12, device = "pdf", path = "/Users/brenna/Dropbox (Partners HealthCare)/MS-chromatin_fine_mapping/results/Fine-mapping/TEAD2_prioritized_Brenna/")
ggsave("common.enrichmentlabeleddots.FullTFT.pdf", p13, device = "pdf", path = "/Users/brenna/Dropbox (Partners HealthCare)/MS-chromatin_fine_mapping/results/Fine-mapping/TEAD2_prioritized_Brenna/")
ggsave("cd4Un.enrichmentlabeleddots.FullTFT.pdf", p14, device = "pdf", path = "/Users/brenna/Dropbox (Partners HealthCare)/MS-chromatin_fine_mapping/results/Fine-mapping/TEAD2_prioritized_Brenna/")
ggsave("BUn.enrichmentlabeleddots.FullTFT.pdf", p15, device = "pdf", path = "/Users/brenna/Dropbox (Partners HealthCare)/MS-chromatin_fine_mapping/results/Fine-mapping/TEAD2_prioritized_Brenna/")

```
