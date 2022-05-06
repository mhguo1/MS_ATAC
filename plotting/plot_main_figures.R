require(RColorBrewer)
require(ggplot2)
require(stringr)
require(locfdr)
require(corrplot)
require(gplots)

#Make Buenrostro plotting parameters
buen.cells<-c("HSC" , "MPP", "LMPP", "CMP",  "CLP", "B", "CD4", "CD8","NK", "pDC", "GMP", "Mono","mDC", "MEP","Ery", "Mega")
pastel<-brewer.pal(4, "Set2")
buen.colors<-c(rep(pastel[1],4), rep(pastel[2], 6), rep(pastel[3], 3), rep(pastel[4], 3))
buen.alphas<-c(seq(0.5, 1, length.out=4), seq(0.5, 1, length.out=6),seq(0.5, 1, length.out=3), seq(0.5, 1, length.out=3))


###Read in Buenrostro Data
#Read in data for single model
buen.single<-do.call(rbind,lapply(list.files("/results/", pattern="single_buenrostro_ATAC", full.names = T),
                                 read.delim))
buen.single<-buen.single[grepl("atac", buen.single$Category, ignore.case=T),]
buen.single$Category<-str_split_fixed(buen.single$Category, "\\_", 3)[,1]

#Read in data for joint model
buen.joint<-do.call(rbind,lapply(list.files("/results/", pattern="joint_buenrostro_ATAC", full.names = T),
                                  read.delim))
buen.joint<-buen.joint[grepl("atac", buen.joint$Category, ignore.case=T),]
buen.joint$Category<-str_split_fixed(buen.joint$Category, "\\_", 3)[,1]

#Read in data for pairwise model
buen.pairwise<-do.call(rbind,
             lapply(list.files("/results/LDSC/Buenrostro/pairwise", pattern="pairwise_buenrostro_ATAC", full.names = T),
                    read.delim))
buen.pairwise$sample<-str_split_fixed(buen.pairwise$Name, "\\_d", 2)[,1]
buen.pairwise$delta<-str_split_fixed(buen.pairwise$Name, "\\_d", 2)[,2]
buen.pairwise<-subset(buen.pairwise, sample%in%buen.cells & delta%in%buen.cells)


#Figure 1B: Single model plotting
dat.plot<-buen.single
dat.plot$Category<-factor(dat.plot$Category, levels=rev(buen.cells))
dat.plot<-dat.plot[order(dat.plot$Category),]
dat.plot$colors<-buen.colors
dat.plot$alphas<-buen.alphas
dat.plot$sizes<-dat.plot$Prop._h2
dat.plot$sizes<-1
ggplot(dat.plot, aes(x=Category, y=-log10(Enrichment_p),  size=sizes))+
  geom_point(col=rev(buen.colors), alpha=rev(buen.alphas), shape=16)+
  coord_flip()+theme_bw()+geom_hline(yintercept = -log10(0.05/nrow(dat.plot)),linetype="dashed")

#Figure 2A: Joint model plotting
dat.plot<-merge(subset(buen.joint, select=c(Category, Coefficient_z.score)), subset(buen.single, select=c(Category, Enrichment_p)), by="Category")
dat.plot$Coefficient_p<-pnorm(-as.numeric(dat.plot$Coefficient_z.score)) #Calculate p values from stratified z-scores
dat.plot$Category<-factor(dat.plot$Category, levels=rev(buen.cells))
dat.plot<-dat.plot[order(dat.plot$Category),]
dat.plot$sizes<-0.7*(-log10(dat.plot$Enrichment_p)) #Calculate sizes of circles based on stratified enrichment p-values
ggplot(dat.plot, aes(x=Category, y=-log10(Coefficient_p), size=sizes)) + 
  geom_point(col=rev(buen.colors), alpha=rev(buen.alphas)) +   
  coord_flip()+theme_bw()+geom_hline(yintercept = -log10(0.05/nrow(dat.plot)),linetype="dashed")


#Figure 2B pairwise model
#Merge pre-conditional p-values
dat.plot<-merge(buen.pairwise, subset(buen.single, select=c(Category, Enrichment_p)), by.x="sample",by.y= "Category", all.x=T,all.y=T)
dat.plot$label<-paste(buen.pairwise$sample, " (", formatC(dat.plot$Enrichment_p, format = "e", digits = 2), ")", sep="")

#Make square matrix with pairwise p-values
corr.dat<-matrix(data=0,nrow=length(buen.cells),ncol=length(buen.cells))
rownames(corr.dat)<-buen.cells
colnames(corr.dat)<-buen.cells

for(i in 1:nrow(corr.dat)){
  for(j in 1:ncol(corr.dat)){
    s1=buen.cells[i]
    s2=buen.cells[j]
    if(i!=j){
     corr.dat[i,j]<-(-log10(subset(dat.plot, sample==s1 & delta==s2)$Coefficient_P_value))
 
    }
    rownames(corr.dat)[i]<-subset(dat.plot, sample==s1)$label[1]
    colnames(corr.dat)[j]<-buen.cells[j]
  }
}
corrplot(corr.dat,p.mat=(10^(-corr.dat)),method = "color", addgrid.col = "black",pch.col = "red",pch.cex=2, type="full", sig.level=0.05/((nrow(corr.dat)-1)^2),insig = "label_sig",is.corr=F, diag=F, cl.lim=c(0,20), tl.col="black")




##Read in Calderon data
cd4.cells<-c("Naive_Teffs_U","Th1_precursors_U", "Th2_precursors_U", "Th17_precursors_U", "Follicular_T_Helper_U","Naive_Tregs_U", "Memory_Tregs_U" )
b.cells<-c("Naive_B_U","Mem_B_U", "Plasmablasts_U")
cd4.cells<-c("Naive_Teffs","Th1_precursors", "Th2_precursors", "Th17_precursors", "Follicular_T_Helper","Naive_Tregs", "Memory_Tregs" )
b.cells<-c("Naive_B","Mem_B", "Plasmablasts")

#Make Calderon plotting parameters
pastel<-brewer.pal(3, "GnBu")
cd4.colors<-c(rep(pastel[2], 5), rep(pastel[3], 2))
cd4.alphas<-c(seq(0.5, 1, length.out=6)[c(1,3:6)],seq(0.5,1,length.out=2))


b.colors<-rep(pastel[3], length(b.cells))
b.alphas<-seq(0.5, 1, length.out=length(b.cells))

cald.single<-do.call(rbind,lapply(list.files("/results/LDSC/Calderon/single", pattern="single_calderon_ATAC", full.names = T),
                                  read.delim))
cald.single<-cald.single[grepl("atac", cald.single$Category, ignore.case=T),]
cald.single$Category<-str_split_fixed(cald.single$Category, "\\_U", 2)[,1]

cald.joint<-do.call(rbind,lapply(list.files("/results", pattern="joint_calderon_ATAC", full.names = T),
                                 read.delim))
cald.joint<-cald.joint[grepl("atac", cald.joint$Category, ignore.case=T),]
cald.joint$Category<-str_split_fixed(cald.joint$Category, "\\_U", 2)[,1]

cald.pairwise<-do.call(rbind,
                       lapply(list.files("/results/", pattern="pairwise_calderon_ATAC", full.names = T),
                              read.delim))

cald.pairwise$sample<-str_split_fixed(cald.pairwise$Name, "\\_d", 2)[,1]
cald.pairwise$delta<-str_split_fixed(cald.pairwise$Name, "\\_d", 2)[,2]
cald.pairwise$sample<-gsub("_U", "", cald.pairwise$sample)
cald.pairwise$delta<-gsub("_U", "", cald.pairwise$delta)




#Figure 3B: Single model plotting
dat.plot<-subset(cald.single, Category%in%cd4.cells)
dat.plot$Category<-factor(dat.plot$Category, levels=rev(cd4.cells))
dat.plot<-dat.plot[order(dat.plot$Category),]
dat.plot$colors<-cd4.colors
dat.plot$alphas<-cd4.alphas
dat.plot$sizes<-dat.plot$Prop._h2
dat.plot$sizes<-1
ggplot(dat.plot, aes(x=Category, y=-log10(Enrichment_p),  size=sizes))+
  geom_point(col=rev(cd4.colors), alpha=rev(cd4.alphas), shape=16)+
  coord_flip()+theme_bw()+geom_hline(yintercept = -log10(0.05/nrow(dat.plot)),linetype="dashed")


#Figure 3c: Joint model plotting
dat.plot<-merge(subset(cald.joint, select=c(Category, Coefficient_z.score)), subset(cald.single, select=c(Category, Enrichment_p)), by="Category")
dat.plot<-subset(dat.plot, Category%in%cd4.cells)
dat.plot$Coefficient_p<-pnorm(-as.numeric(dat.plot$Coefficient_z.score)) #Calculate p values from stratified z-scores
dat.plot$Category<-factor(dat.plot$Category, levels=rev(cd4.cells))
dat.plot<-dat.plot[order(dat.plot$Category),]
dat.plot$sizes<-0.7*(-log10(dat.plot$Enrichment_p)) #Calculate sizes of circles based on stratified enrichment p-values
ggplot(dat.plot, aes(x=Category, y=-log10(Coefficient_p), size=sizes)) + 
  geom_point(col=rev(cd4.colors), alpha=rev(cd4.alphas)) +   
  coord_flip()+theme_bw()+geom_hline(yintercept = -log10(0.05/nrow(dat.plot)),linetype="dashed")


#Figure 4B: Single model plotting
dat.plot<-subset(cald.single, Category%in%b.cells)
dat.plot$Category<-factor(dat.plot$Category, levels=rev(b.cells))
dat.plot<-dat.plot[order(dat.plot$Category),]
dat.plot$colors<-b.colors
dat.plot$alphas<-b.alphas
dat.plot$sizes<-1ggplot(dat.plot, aes(x=Category, y=-log10(Enrichment_p),  size=sizes))+
  geom_point(col=rev(b.colors), alpha=rev(b.alphas), shape=16)+
  coord_flip()+theme_bw()+geom_hline(yintercept = -log10(0.05/nrow(dat.plot)),linetype="dashed")


#Figure 4c: Joint model plotting
dat.plot<-merge(subset(cald.joint, select=c(Category, Coefficient_z.score)), subset(cald.single, select=c(Category, Enrichment_p)), by="Category")
dat.plot<-subset(dat.plot, Category%in%b.cells)
dat.plot$Coefficient_p<-pnorm(-as.numeric(dat.plot$Coefficient_z.score)) #Calculate p values from stratified z-scores
dat.plot$Category<-factor(dat.plot$Category, levels=rev(b.cells))
dat.plot<-dat.plot[order(dat.plot$Category),]
dat.plot$sizes<-0.7*(-log10(dat.plot$Enrichment_p)) #Calculate sizes of circles based on stratified enrichment p-values
ggplot(dat.plot, aes(x=Category, y=-log10(Coefficient_p), size=sizes)) + 
  geom_point(col=rev(b.colors), alpha=rev(b.alphas)) +   
  coord_flip()+theme_bw()+geom_hline(yintercept = -log10(0.05/nrow(dat.plot)),linetype="dashed")




###Verily Data###
#Make Verily plotting parameters
verily.cells<-c("traB", "cMBc", "T4nv", "T4cm", "T4em", "T4ra")
pastel<-brewer.pal(2, "Set2")
verily.colors<-c(rep(pastel[1],2), rep(pastel[2], 4))
verily.alphas<-c(seq(0.5, 1, length.out=2), seq(0.5, 1, length.out=4))


#Verily untreated, single model
#Read in data for single model
verily.single<-do.call(rbind,lapply(list.files("/results/Verily/untreated/", pattern="single_verily_ATAC", full.names = T),
                                  read.delim))
verily.single<-verily.single[grepl("untreated", verily.single$Category, ignore.case=T),]
verily.single$Category<-str_split_fixed(verily.single$Category, "\\_", 3)[,1]

#Figure 5a: Single model plotting
dat.plot<-verily.single
dat.plot$Category<-factor(dat.plot$Category, levels=rev(verily.cells))
dat.plot<-dat.plot[order(dat.plot$Category),]
dat.plot$colors<-verily.colors
dat.plot$alphas<-verily.alphas
dat.plot$sizes<-1
ggplot(dat.plot, aes(x=Category, y=-log10(Enrichment_p),  size=sizes))+
  geom_point(col=rev(verily.colors), alpha=rev(verily.alphas), shape=16)+
  coord_flip(ylim=c(0,4))+theme_bw()+geom_hline(yintercept = -log10(0.05/nrow(dat.plot)),linetype="dashed")


#Figure 5b: Joint model, CD4 T cells verily
#Verily untreated, joint model with CD4
#Read in data for joint model
verily.cd4.joint<-read.delim("/results/Verily/untreated/CD4/stratified.results", header=T, stringsAsFactors = F, sep="\t")
verily.cd4.joint<-verily.cd4.joint[grepl("untreated", verily.cd4.joint$Category, ignore.case=T),]
verily.cd4.joint$Category<-str_split_fixed(verily.cd4.joint$Category, "\\_", 3)[,1]

verily.cd4.cells<-verily.cells[3:6]
verily.cd4.alphas<-verily.alphas[3:6]
verily.cd4.colors<-verily.colors[3:6]

dat.plot<-merge(subset(verily.cd4.joint, select=c(Category, Coefficient_z.score)), subset(verily.single, select=c(Category, Enrichment_p)), by="Category")
dat.plot$Coefficient_p<-pnorm(-as.numeric(dat.plot$Coefficient_z.score)) #Calculate p values from stratified z-scores
dat.plot$Category<-factor(dat.plot$Category, levels=rev(verily.cd4.cells))
dat.plot<-dat.plot[order(dat.plot$Category),]
dat.plot$sizes<-0.7*(-log10(dat.plot$Enrichment_p)) #Calculate sizes of circles based on stratified enrichment p-values
ggplot(dat.plot, aes(x=Category, y=-log10(Coefficient_p), size=sizes)) + 
  geom_point(col=rev(verily.cd4.colors), alpha=rev(verily.cd4.alphas)) +   
  coord_flip(ylim=c(0,3))+theme_bw()+geom_hline(yintercept = -log10(0.05/nrow(dat.plot)),linetype="dashed")



#Figure 5C: Joint model, B cells verily
#Verily untreated, joint model with B
#Read in data for joint model
verily.b.joint<-read.delim("/results/Verily/untreated/B/stratified.results", header=T, stringsAsFactors = F, sep="\t")
verily.b.joint<-verily.b.joint[grepl("untreated", verily.b.joint$Category, ignore.case=T),]
verily.b.joint$Category<-str_split_fixed(verily.b.joint$Category, "\\_", 3)[,1]

verily.b.cells<-verily.cells[1:2]
verily.b.alphas<-verily.alphas[1:2]
verily.b.colors<-verily.colors[1:2]

dat.plot<-merge(subset(verily.b.joint, select=c(Category, Coefficient_z.score)), subset(verily.single, select=c(Category, Enrichment_p)), by="Category")
dat.plot$Coefficient_p<-pnorm(-as.numeric(dat.plot$Coefficient_z.score)) #Calculate p values from stratified z-scores
dat.plot$Category<-factor(dat.plot$Category, levels=rev(verily.b.cells))
dat.plot<-dat.plot[order(dat.plot$Category),]
dat.plot$sizes<-0.7*(-log10(dat.plot$Enrichment_p)) #Calculate sizes of circles based on stratified enrichment p-values
ggplot(dat.plot, aes(x=Category, y=-log10(Coefficient_p), size=sizes)) + 
  geom_point(col=rev(verily.b.colors), alpha=rev(verily.b.alphas)) +   
  coord_flip(ylim=c(0,3))+theme_bw()+geom_hline(yintercept = -log10(0.05/nrow(dat.plot)),linetype="dashed")



#Figure 5D #Treatment enrichments
verily.natalizumab.single<-do.call(rbind,lapply(list.files("/results/Verily/natalizumab/", pattern=".results", full.names = T),
                                    read.delim))
verily.natalizumab.single$treatment<-"natalizumab"
verily.glatiramer.single<-do.call(rbind,lapply(list.files("/results/Verily/glatiramer/", pattern=".results", full.names = T),
                                    read.delim))
verily.glatiramer.single$treatment<-"glatiramer"
verily.interferon.single<-do.call(rbind,lapply(list.files("/results/Verily/interferon/", pattern=".results", full.names = T),
                                     read.delim))
verily.interferon.single$treatment<-"interferon"
verily.treated.single<-rbind(verily.natalizumab.single, verily.glatiramer.single, verily.interferon.single)

verily.treated.single<-verily.treated.single[grepl("_1", verily.treated.single$Category, ignore.case=T),]
verily.treated.single$Category<-str_split_fixed(verily.treated.single$Category, "\\_", 3)[,1]

dat.plot<-subset(verily.treated.single, treatment=="interferon")
dat.plot$Category<-factor(dat.plot$Category, levels=rev(verily.cells))
dat.plot<-dat.plot[order(dat.plot$Category),]
dat.plot$colors<-verily.colors
dat.plot$alphas<-verily.alphas
dat.plot$sizes<-1
ggplot(dat.plot, aes(x=Category, y=-log10(Enrichment_p),  size=sizes))+
  geom_point(col=rev(verily.colors), alpha=rev(verily.alphas), shape=16)+
  coord_flip(ylim=c(0,4))+theme_bw()+geom_hline(yintercept = -log10(0.05/nrow(dat.plot)),linetype="dashed")


