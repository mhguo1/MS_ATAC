## List creation

# cd to the directory with the CD4 and B gene lists
#cd /Users/np911/Dropbox (Partners HealthCare)/Projects/MS-chromatin_fine_mapping/results/Fine-mapping 

# CD4 
wc -l CD4_prioritized_genes.txt # 364
wc -l B_prioritized_genes.txt # 261

# find the common gene list
cat CD4_prioritized_genes.txt B_prioritized_genes.txt | awk 'x[$0]++' | sort > common_prioritized_genes.txt
wc -l common_prioritized_genes.txt # 178

## create the unique gene lists
R --no-restore --no-save 
cd4 <- read.table("CD4_prioritized_genes.txt", header=FALSE, stringsAsFactors=FALSE)
b <- read.table("B_prioritized_genes.txt", header=FALSE, stringsAsFactors=FALSE)
common <- read.table("common_prioritized_genes.txt", header=FALSE, stringsAsFactors=FALSE)

## CD4 unique
sum(!(cd4$V1 %in% common$V1))
# 186
cd4_unique <- as.data.frame(cd4$V1[!(cd4$V1 %in% common$V1)])
write.table(cd4_unique, file="cd4_unique_prioritized_genes.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

## B unique
sum(!(b$V1 %in% common$V1))
# 83
b_unique <- as.data.frame(b$V1[!(b$V1 %in% common$V1)])
write.table(b_unique, file="b_unique_prioritized_genes.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)
q()
