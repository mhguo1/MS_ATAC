#!/usr/bin/sh


#Step 1: Download data from NCBI GEO GSE118189
/software/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump -I --split-files ${SRR}
mv ${SRR}_1.fastq ${sample}_1.fastq
mv ${SRR}_2.fastq ${sample}_2.fastq

#Step 2: Trim adapters
/software/cutadapt -m 20 -O 5 -a CTGTCTCTTATA -A CTGTCTCTTATA -o ${sample}.trim.1.fastq -p ${sample}.trim.2.fastq ${sample}_1.fastq ${sample}_2.fastq | gzip -c > ${sample}.cutadapt.log.gz
gzip ${sample}.trim.1.fastq
gzip ${sample}.trim.2.fastq

#Step 3: Align using HISAT2
hisat2 -x /hisat2_index/grch37_snp/genome_snp -1 ${sample}.trim.1.fastq.gz -2 ${sample}.trim.2.fastq.gz | samtools view -hb - > ${sample}.hisat2.aln.bam
samtools sort -o ${sample}.sort.bam -T ${sample}.tmp ${sample}.hisat2.aln.bam
samtools index ${sample}.sort.bam

#Step 4: Filter bam files
samtools view -hS -q 30 -F 1804 -f 2 ${sample}.sort.bam | awk '$3!="MT"' | samtools view -hb - > ${sample}.sort.filter.bam
samtools index ${sample}.sort.filter.bam

java -Xmx12G -jar /software//bin/picard-private.jar MarkDuplicates INPUT=${sample}.sort.filter.bam OUTPUT=${sample}.dedup.bam METRICS_FILE=${sample}.picard_metrics.txt REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT
samtools index ${sample}.dedup.bam

#Step 5: Samples from same cell type were then merged into a single bam
samtools merge ${prefix}.merge.dedup.bam ${prefix}_*.dedup.bam
samtools index ${prefix}.merge.dedup.bam

#Step 5: Call peaks on merged bam files
macs2 callpeak -t ${sample}.dedup.bam -f BAM -n ${sample}.MACS -g hs --keep-dup all -q 0.1 --nomodel --nolambda
