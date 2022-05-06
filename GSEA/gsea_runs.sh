## Core to replicate gsea enrichment analyses (all these are local on Niko's laptop)
## Update: all were re-run Feb 12th 2021 with updated final lists. The previous results were moved to "archive" subdirectory.


cd /Users/np911/Dropbox \(Partners HealthCare\)/Projects/MS-chromatin_fine_mapping/results/Fine-mapping

## -- CD4+ T cells -- ##

mkdir -p gsea/CD4_T
for i in h c1 cgp cp mirdb gtrd c5 c6 c7
do
 #python2.7 /Users/np911/Documents/Projects/InProgress/Cmap/gsea/gsea-v2.py -i CD4_T_cells_genes.txt  -q 1 "$i"  >  gsea/CD4_T/$i.enrichment
 python2.7 /Users/np911/Documents/Projects/InProgress/Cmap/gsea/gsea-v2.py -i CD4_prioritized_genes.txt  -q 1 "$i"  >  gsea/CD4_T/$i.enrichment
 mv gseaPlot1.pdf  gsea/CD4_T/gseaPlot1.$i.pdf
 mv gseaPlot2.pdf  gsea/CD4_T/gseaPlot2.$i.pdf
done


## -- B cells -- ##

mkdir -p gsea/B
for i in h c1 cgp cp mirdb gtrd c5 c6 c7
do
 #python2.7 /Users/np911/Documents/Projects/InProgress/Cmap/gsea/gsea-v2.py -i B_cells_genes.txt  -q 1 "$i"  >  gsea/B/$i.enrichment
 python2.7 /Users/np911/Documents/Projects/InProgress/Cmap/gsea/gsea-v2.py -i B_prioritized_genes.txt -q 1 "$i"  >  gsea/B/$i.enrichment
 mv gseaPlot1.pdf  gsea/B/gseaPlot1.$i.pdf
 mv gseaPlot2.pdf  gsea/B/gseaPlot2.$i.pdf
done

## -- Common genes -- ##

mkdir -p gsea/CD4_B_common
for i in h c1 cgp cp mirdb gtrd c5 c6 c7
do
 #python2.7 /Users/np911/Documents/Projects/InProgress/Cmap/gsea/gsea-v2.py -i common_genes.txt  -q 1 "$i"  >  gsea/CD4_B_common/$i.enrichment
 python2.7 /Users/np911/Documents/Projects/InProgress/Cmap/gsea/gsea-v2.py -i common_prioritized_genes.txt  -q 1 "$i"  >  gsea/CD4_B_common/$i.enrichment
 mv gseaPlot1.pdf  gsea/CD4_B_common/gseaPlot1.$i.pdf
 mv gseaPlot2.pdf  gsea/CD4_B_common/gseaPlot2.$i.pdf
done

## -- Unique CD4 T genes -- ##

mkdir -p gsea/CD4_unique
for i in h c1 cgp cp mirdb gtrd c5 c6 c7
do
 #python2.7 /Users/np911/Documents/Projects/InProgress/Cmap/gsea/gsea-v2.py -i unique_CD4_genes.txt  -q 1 "$i"  >  gsea/CD4_unique/$i.enrichment
 python2.7 /Users/np911/Documents/Projects/InProgress/Cmap/gsea/gsea-v2.py -i cd4_unique_prioritized_genes.txt  -q 1 "$i"  >  gsea/CD4_unique/$i.enrichment
 mv gseaPlot1.pdf  gsea/CD4_unique/gseaPlot1.$i.pdf
 mv gseaPlot2.pdf  gsea/CD4_unique/gseaPlot2.$i.pdf
done



## -- Unique B T genes -- ##

mkdir -p gsea/B_unique
for i in h c1 cgp cp mirdb gtrd c5 c6 c7
do
 #python2.7 /Users/np911/Documents/Projects/InProgress/Cmap/gsea/gsea-v2.py -i unique_B_genes.txt  -q 1 "$i"  >  gsea/B_unique/$i.enrichment
 python2.7 /Users/np911/Documents/Projects/InProgress/Cmap/gsea/gsea-v2.py -i b_unique_prioritized_genes.txt  -q 1 "$i"  >  gsea/B_unique/$i.enrichment
 mv gseaPlot1.pdf  gsea/B_unique/gseaPlot1.$i.pdf
 mv gseaPlot2.pdf  gsea/B_unique/gseaPlot2.$i.pdf
done
