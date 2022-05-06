library("readxl")
###Read in Table S18 from Science paper
setwd("/Dropbox/MS-chromatin_fine_mapping/results/Fine-mapping/")

#table<-read.delim("MS_Science_Table_S18.txt", header=T, stringsAsFactors = F, sep="\t")
table<-data.frame(read_xlsx("NIHMS1580882-supplement-Supplementary_Tables_11-20.xlsx", sheet=8,skip=4), stringsAsFactors = F)
exonic_genes<-unique(unlist(strsplit(paste(table$Exonic.genes[!is.na(table$Exonic.genes)], collapse="|"), "\\|")))
eqtl_genes<-unique(unlist(strsplit(paste(table$eQTL.genes[!is.na(table$eQTL.genes)], collapse="|"), "\\|")))
network_genes<-unique(unlist(strsplit(paste(table$Regulatory.network[!is.na(table$Regulatory.network)], collapse="|"), "\\|")))
depict_genes<-unique(unlist(strsplit(paste(table$Depict[!is.na(table$Depict)], collapse="|"), "\\|")))

ms2019_gene_list<-unique(c(exonic_genes, eqtl_genes, network_genes, depict_genes))
ms2019_gene_list<-ms2019_gene_list[!is.na(ms2019_gene_list)]
#length(unique(c(exonic_genes, eqtl_genes, network_genes, depict_genes))) #I get 560 genes


##Generate master table with prioritized genes
setwd("/Dropbox/MS-chromatin_fine_mapping/results/Fine-mapping/")

dat<-read.delim("MS.PICS.ld0.2.ATAC.gencode_v19.PCHiC.motifbreakr.eqtl.cs.bed", header=T, stringsAsFactors = F, sep="\t")

#Add in PCHiC data for B and CD4
dat$PCHiC_B_genes<-NA
dat$PCHiC_CD4_genes<-NA

for(i in 1:nrow(dat)){
  #Parse out CD4 PCHiC looped genes
  if(grepl("CD4", dat[i,]$pchic_cells) | grepl("B", dat[i,]$pchic_cells)){
    temp_cells<-unlist(strsplit(dat[i,]$pchic_cells, "\\|"))
    genes_cd4<-vector()
    for(c in 1:length(temp_cells)){
      if(grepl("CD4", temp_cells[c])){
        genes_cd4<-c(genes_cd4, unlist(strsplit( unlist(strsplit(dat[i,]$pchic_genes, "\\|"))[c],"\\;")))
       }
    }
    if(length(genes_cd4)>0){
      dat[i,]$PCHiC_CD4_genes<-paste(genes_cd4, collapse=",")
    }
    
    genes_b<-vector()
    for(c in 1:length(temp_cells)){
      if(grepl("B", temp_cells[c])){
        genes_b<-c(genes_b, unlist(strsplit( unlist(strsplit(dat[i,]$pchic_genes, "\\|"))[c],"\\;")))
      }
    }
    if(length(genes_b)>0){
      dat[i,]$PCHiC_B_genes<-paste(genes_b, collapse=",")
    }
  }
}


#Parse out PCHiC targets for B cells
#Create empty table, where each locus is its own row
dat.summary<-data.frame(lead_snp=unique(dat$lead_snp), ATAC_B=0, ATAC_CD4=0, ATAC_CD4_subset=0, ATAC_Th17=0, ATAC_B_subset=0 , ATAC_MemB=0, 
                        PCHiC_CD4=NA, PCHiC_B=NA, CD4_genes=NA, B_genes=NA, CS_snps=0, stringsAsFactors = F )
for(i in 1:nrow(dat.summary)){
  
  #temp.dat is a temporary df for each locus
  temp.dat<-subset(dat, lead_snp==dat.summary[i,]$lead_snp & CS_95==1 & pics>0.01)
 # dat.summary[i,]$locus<-dat[dat$lead_snp==dat.summary[i,]$lead_snp,]$locus[1]
  
  if(nrow(temp.dat)>0){
    #Parse out information for that locus
    dat.summary[i,]$ATAC_B<-max(temp.dat$B)
    dat.summary[i,]$ATAC_CD4<-max(temp.dat$CD4)
    dat.summary[i,]$ATAC_Th17<-max(temp.dat$Th17_precursors_U)
    dat.summary[i,]$ATAC_CD4_subset<-max(temp.dat$Memory_Teffs_U, temp.dat$Memory_Tregs_U, temp.dat$Naive_Teffs_U, temp.dat$Naive_Tregs_U, temp.dat$Th17_precursors_U, temp.dat$Follicular_T_Helper_U, temp.dat$Th1_precursors_U, temp.dat$Th2_precursors_U)
    dat.summary[i,]$ATAC_MemB<-max(temp.dat$Mem_B_U)
    dat.summary[i,]$ATAC_B_subset<-max(temp.dat$Mem_B_U, temp.dat$Naive_B_U, temp.dat$Plasmablasts_U)
    dat.summary[i,]$CS_snps<-nrow(temp.dat)
    
  
    if(length(temp.dat$PCHiC_B_genes[!is.na(temp.dat$PCHiC_B_genes)])>0){
      dat.summary[i,]$PCHiC_B<-paste(unique(unlist(strsplit(paste(temp.dat$PCHiC_B_genes[!is.na(temp.dat$PCHiC_B_genes)], collapse=","), "\\,"))), collapse=",")
    }
    if(length(temp.dat$PCHiC_CD4_genes[!is.na(temp.dat$PCHiC_CD4_genes)])>0){
      dat.summary[i,]$PCHiC_CD4<-paste(unique(unlist(strsplit(paste(temp.dat$PCHiC_CD4_genes[!is.na(temp.dat$PCHiC_CD4_genes)], collapse=","), "\\,"))), collapse=",")
    }
    b_genes<-unique(unlist(strsplit(paste(subset(temp.dat, !is.na(PCHiC_B_genes) & ( B==1 | Mem_B_U==1 | Naive_B_U==1 | Plasmablasts_U==1 ))$PCHiC_B_genes, collapse=","), "\\,")))
    cd4_genes<-unique(unlist(strsplit(paste(subset(temp.dat, !is.na(PCHiC_CD4_genes) & ( CD4==1 | Th17_precursors_U==1 | Memory_Teffs_U==1 | Memory_Tregs_U==1 | Memory_Tregs_U==1 | Naive_Teffs_U==1 | Naive_Tregs_U==1 | Follicular_T_Helper_U==1 | Th1_precursors_U==1 | Th2_precursors_U==1 ))$PCHiC_CD4_genes, collapse=","), "\\,")))
    if(length(cd4_genes)>0){
      dat.summary[i,]$CD4_genes<-paste0(cd4_genes, collapse=",")
    }
    if(length(b_genes)>0){
      dat.summary[i,]$B_genes<-paste0(b_genes, collapse=",")
    }
   

  }
}
dat.summary$lead_chr<-gsub("chr", "", str_split_fixed(dat.summary$lead_snp, "\\:", 4)[,1])
dat.summary$lead_pos<-str_split_fixed(dat.summary$lead_snp, "\\:", 4)[,2]
write.table(dat.summary, "MS_locus_prioritization.txt", row.names = F, col.names = T, sep="\t", quote=F)

cd4_gene_list<-unique(unlist(strsplit(dat.summary$CD4_genes[!is.na(dat.summary$CD4_genes)], "\\,")))
write.table(cd4_gene_list, "CD4_prioritized_genes.txt", row.names = F, col.names = F, sep="\t", quote=F)
b_gene_list<-unique(unlist(strsplit(dat.summary$B_genes[!is.na(dat.summary$B_genes)], "\\,")))
write.table(b_gene_list, "B_prioritized_genes.txt", row.names = F, col.names = F, sep="\t", quote=F)

length(unique(ms2019_gene_list[!is.na(ms2019_gene_list)])) #560 genes
length(intersect(cd4_gene_list, ms2019_gene_list))#121 genes shared between our paper and the 2019 Science paper for CD4
length(intersect(b_gene_list, ms2019_gene_list)) #67 genes shared between our paper and the 2019 Science paper for CD4

