library("exomePeak")
setwd("results/exomepeak")
INPUT=c("results/bam/293T-IPIN.filter.sorted.bam")
IP=c("results/bam/293T-IP1.filter.sorted.bam")
GENE_ANNO_GTF="hg38/gencode.v31.annotation.gtf"
result=exomepeak(GENE_ANNO_GTF=GENE_ANNO_GTF, MINIMAL_MAPQ = 30, IP_BAM=IP, INPUT_BAM=INPUT,  EXPERIMENT_NAME="IP1", WINDOW_WIDTH=200)
