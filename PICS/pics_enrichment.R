#Enrichment of PICS versus ATAC
dat.cs<-read.delim("MS.PICS.ld0.2.ATAC.gencode_v19.PCHiC.motifbreakr.eqtl.cs.bed", header=T, stringsAsFactors = F, sep="\t")
#cell<-c("B","CD4","CD8","CLP","CMP","Ery","GMP","HSC","LMPP","mDC","MEP","MPP","Mega","Mono","NK","pDC")    
#cell<-c("Central_memory_CD8pos_T_U","Effector_memory_CD8pos_T_U","Gamma_delta_T_U","Immature_NK_U","Mature_NK_U","Mem_B_U","Memory_Teffs_U","Memory_Tregs_U","Monocytes_U","Myeloid_DCs_U","Naive_B_U","Naive_CD8_T_U","Naive_Teffs_U","Naive_Tregs_U","pDCs_U","Plasmablasts_U","Th17_precursors_U","Follicular_T_Helper_U","Memory_NK_U","Th1_precursors_U","Th2_precursors_U")
cell<-c("Mem_B_U","Memory_Tregs_U","Naive_B_U","Naive_Teffs_U","Naive_Tregs_U","Plasmablasts_U","Th17_precursors_U","Follicular_T_Helper_U","Th1_precursors_U","Th2_precursors_U")

dat.results<-data.frame(cell_type=cell, fold_change=0, p_value=0, stringsAsFactors = F)

for(i in 1:nrow(dat.results)){
  dat.cs$temp_col<-dat.cs[,which(names(dat.cs)==dat.results[i,]$cell_type)]
  a<-nrow(subset(dat.cs, (CS_95==1 & pics>0.01) & temp_col==1))
  b<-nrow(subset(dat.cs, (CS_95==1 & pics>0.01) & temp_col==0))
  c<-nrow(subset(dat.cs, !(CS_95==1 & pics>0.01) & temp_col==1))
  d<-nrow(subset(dat.cs, !(CS_95==1 & pics>0.01) & temp_col==0))
  
  
  dat.results[i,]$fold_change<-(a/(a+b))/(c/(c+d))
  dat.results[i,]$p_value<-fisher.test(rbind(c(a,b), c(c, d)), alternative="greater")$p.value
}
dat.results[order(dat.results$p_value),]
write.table(dat.results, "MS_ATAC_PICS_enrichment.txt", row.names = F, col.names = T, quote=F, sep="\t")


#Pairwise PICS ATAC Enrichment
dat.cs<-read.delim("MS.PICS.ld0.2.ATAC.gencode_v19.PCHiC.motifbreakr.eqtl.cs.bed", header=T, stringsAsFactors = F, sep="\t")
cell<-c("Mem_B_U","Memory_Tregs_U","Naive_B_U","Naive_Teffs_U","Naive_Tregs_U","Plasmablasts_U","Th17_precursors_U","Follicular_T_Helper_U","Th1_precursors_U","Th2_precursors_U")

#Make square matrix with pairwise p-values
corr.dat<-matrix(data=0,nrow=length(cell),ncol=length(cell))
rownames(corr.dat)<-cell
colnames(corr.dat)<-cell

for(i in 1:nrow(corr.dat)){
  for(j in 1:ncol(corr.dat)){
  s1<-cell[i]
  s2<-cell[j]
  dat.cs$temp_col<-dat.cs[,which(names(dat.cs)==s1) ]-dat.cs[,which(names(dat.cs)==s2)]
  a<-nrow(subset(dat.cs, (CS_95==1 & pics>0.01) & temp_col==1))
  b<-nrow(subset(dat.cs, (CS_95==1 & pics>0.01) & temp_col!=1))
  c<-nrow(subset(dat.cs, !(CS_95==1 & pics>0.01) & temp_col==1))
  d<-nrow(subset(dat.cs, !(CS_95==1 & pics>0.01) & temp_col!=1))
  
  corr.dat[i,j]<-fisher.test(rbind(c(a,b), c(c, d)), alternative="greater")$p.value
  }
}
corr.dat<-data.frame(corr.dat)
write.table(corr.dat, "MS_ATAC_PICS_pairwise_enrichment.txt", row.names = T, col.names = T, quote=F, sep="\t")




#Enrichment of PICS versus PCHiC
dat.cs<-read.delim("MS.PICS.ld0.2.ATAC.gencode_v19.PCHiC.motifbreakr.eqtl.cs.bed", header=T, stringsAsFactors = F, sep="\t")
dat.cs<-subset(dat.cs, select=c(locus,pics, CS_95, pchic_cells))
dat.cs$pchic_cells<-gsub("\\|", "\\,",dat.cs$pchic_cells)
dat.cs$pchic_cells<-paste0(",", dat.cs$pchic_cells, ",")
#cell<-c("nB", "tB", "naCD4", "aCD4", "tCD4", "Ery", "FoeT", "EP", "Mac0", "Mac1", "Mac2", "Mon", "MK", "nCD8", "tCD8", "Neu")
cell<-"Mon"
dat.results<-data.frame(cell_type=cell, a=0, b=0, c=0, d=0, fold_change=0, p_value=0, stringsAsFactors = F)

for(i in 1:nrow(dat.results)){
  dat.cs$temp_col<-ifelse(grepl(paste0(",", dat.results[i,]$cell_type, ","), dat.cs$pchic_cells), 1, 0)
  a<-nrow(subset(dat.cs, (CS_95==1 & pics>0.01 ) & temp_col==1 ))
  b<-nrow(subset(dat.cs, (CS_95==1 & pics>0.01) & temp_col==0))
  c<-nrow(subset(dat.cs, !(CS_95==1 & pics>0.01) & temp_col==1))
  d<-nrow(subset(dat.cs, !(CS_95==1 & pics>0.01) & temp_col==0))
  
  dat.results[i,]$a<-a
  dat.results[i,]$b<-b
  dat.results[i,]$c<-c
  dat.results[i,]$d<-d
  dat.results[i,]$fold_change<-(a/(a+b))/(c/(c+d))
  dat.results[i,]$p_value<-fisher.test(rbind(c(a,b), c(c, d)), alternative="greater")$p.value
}
dat.results[order(dat.results$p_value),]
write.table(dat.results, "MS_PCHiC_PICS_enrichment.txt", row.names = F, col.names = T, quote=F, sep="\t")



#Enrichment of PICS versus PCHiC + ATAC-seq
dat.cs<-read.delim("MS.PICS.ld0.2.ATAC.gencode_v19.PCHiC.motifbreakr.eqtl.cs.bed", header=T, stringsAsFactors = F, sep="\t")
dat.cs$pchic_cells<-gsub("\\|", "\\,",dat.cs$pchic_cells)
dat.cs$pchic_cells<-paste0(",", dat.cs$pchic_cells, ",")
cell<-c("nB", "tB", "naCD4", "aCD4", "tCD4", "Ery", "FoeT", "EP", "Mac0", "Mac1", "Mac2", "Mon", "MK", "nCD8", "tCD8", "Neu")
atac<-c("B", "B", "CD4", "CD4", "CD4", "MEP", "B", "MEP", "Mono", "Mono","Mono", "Mono", "Mega", "CD8", "CD8", "GMP" )
dat.results<-data.frame(pchic=cell, atac=atac,a=0, b=0, c=0, d=0, fold_change=0, p_value=0, stringsAsFactors = F)

for(i in 1:nrow(dat.results)){
  dat.cs$temp_col_pchic<-ifelse(grepl(paste0(",", dat.results[i,]$pchic, ","), dat.cs$pchic_cells), 1, 0)
  dat.cs$temp_col_atac<-dat.cs[,which(names(dat.cs)==dat.results[i,]$atac)]
  
  a<-nrow(subset(dat.cs, (CS_95==1 & pics>0.01 ) & (temp_col_pchic==1 & temp_col_atac==1)))
  b<-nrow(subset(dat.cs, (CS_95==1 & pics>0.01) &  !(temp_col_pchic==1 & temp_col_atac==1)))
  c<-nrow(subset(dat.cs, !(CS_95==1 & pics>0.01) &  (temp_col_pchic==1 & temp_col_atac==1)))
  d<-nrow(subset(dat.cs, !(CS_95==1 & pics>0.01) &  !(temp_col_pchic==1 & temp_col_atac==1)))
  
  dat.results[i,]$a<-a
  dat.results[i,]$b<-b
  dat.results[i,]$c<-c
  dat.results[i,]$d<-d
  dat.results[i,]$fold_change<-(a/(a+b))/(c/(c+d))
  dat.results[i,]$p_value<-fisher.test(rbind(c(a,b), c(c, d)), alternative="greater")$p.value
}
dat.results[order(dat.results$p_value),]
dat.results<-subset(dat.results, select=-c(a,b,c,d))
write.table(dat.results, "MS_PCHiC_ATAC_PICS_enrichment.txt", row.names = F, col.names = T, quote=F, sep="\t")





#Enrichment of PICS versus eQTL
dat.cs<-read.delim("MS.PICS.ld0.2.ATAC.gencode_v19.PCHiC.motifbreakr.eqtl.cs.bed", header=T, stringsAsFactors = F, sep="\t")
dat.cs<-subset(dat.cs, select=c(pics, CS_95, eqtl_cells))
#dat.cs$pchic_cells<-gsub("\\|", "\\,",dat.cs$eqtl_cells)
dat.cs$eqtl_cells<-paste0(",", dat.cs$eqtl_cells, ",")
cell<-c("B_CELL_NAIVE","MONOCYTES","M2","NK","TREG_MEM","CD4_NAIVE","CD4_STIM","TREG_NAIVE","TFH","TH1","THSTAR","TH17","TH2","CD8_NAIVE","CD8_STIM")
dat.results<-data.frame(cell_type=cell, fold_change=0, p_value=0, stringsAsFactors = F)

for(i in 1:nrow(dat.results)){
  dat.cs$temp_col<-ifelse(grepl(paste0(",", dat.results[i,]$cell_type, ","), dat.cs$eqtl_cells), 1, 0)
  a<-nrow(subset(dat.cs, (CS_95==1 & pics>0.01) & temp_col==1))
  b<-nrow(subset(dat.cs, (CS_95==1 & pics>0.01) & temp_col==0))
  c<-nrow(subset(dat.cs, !(CS_95==1 & pics>0.01) & temp_col==1))
  d<-nrow(subset(dat.cs, !(CS_95==1 & pics>0.01) & temp_col==0))
  
  dat.results[i,]$fold_change<-(a/(a+b))/(c/(c+d))
  dat.results[i,]$p_value<-fisher.test(rbind(c(a,b), c(c, d)), alternative="greater")$p.value
}
dat.results[order(dat.results$p_value),]



#Pairwise PICS eQTL Enrichment
dat.cs<-read.delim("MS.PICS.ld0.2.ATAC.gencode_v19.PCHiC.motifbreakr.eqtl.cs.bed", header=T, stringsAsFactors = F, sep="\t")
dat.cs<-subset(dat.cs, select=c(pics, CS_95, eqtl_cells))
#dat.cs$pchic_cells<-gsub("\\|", "\\,",dat.cs$eqtl_cells)
dat.cs$eqtl_cells<-paste0(",", dat.cs$eqtl_cells, ",")
cell<-c("B_CELL_NAIVE","MONOCYTES","M2","NK","TREG_MEM","CD4_NAIVE","CD4_STIM","TREG_NAIVE","TFH","TH1","THSTAR","TH17","TH2","CD8_NAIVE","CD8_STIM")

#Make square matrix with pairwise p-values
corr.dat<-matrix(data=0,nrow=length(cell),ncol=length(cell))
rownames(corr.dat)<-cell
colnames(corr.dat)<-cell

for(i in 1:nrow(corr.dat)){
  for(j in 1:ncol(corr.dat)){
    s1<-cell[i]
    s2<-cell[j]
    dat.cs$temp_col_s1<-ifelse(grepl(paste0(",", s1, ","), dat.cs$eqtl_cells), 1, 0)
    dat.cs$temp_col_s2<-ifelse(grepl(paste0(",", s2, ","), dat.cs$eqtl_cells), 1, 0)
    dat.cs$temp_col<-dat.cs$temp_col_s1-dat.cs$temp_col_s2
    
    a<-nrow(subset(dat.cs, (CS_95==1 & pics>0.01) & temp_col==1))
    b<-nrow(subset(dat.cs, (CS_95==1 & pics>0.01) & temp_col!=1))
    c<-nrow(subset(dat.cs, !(CS_95==1 & pics>0.01) & temp_col==1))
    d<-nrow(subset(dat.cs, !(CS_95==1 & pics>0.01) & temp_col!=1))
    
    corr.dat[i,j]<-fisher.test(rbind(c(a,b), c(c, d)), alternative="greater")$p.value
  }
}
corr.dat<-data.frame(corr.dat)
write.table(corr.dat, "MS_eQTL_PICS_pairwise_enrichment.txt", row.names = T, col.names = T, quote=F, sep="\t")

