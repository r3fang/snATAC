## Introduction
snATAC is a Ren-lab in-house analysis pipeline for single-nucleus ATAC-seq (snATAC-seq).

## Requirements

- softwares ([samtools 1.2](http://www.htslib.org/doc/samtools.html), [Python 2.7](https://www.python.org/download/releases/2.7/), [bedtools](http://bedtools.readthedocs.io/en/latest/), [MACS2](https://github.com/taoliu/MACS), [bwa](https://github.com/lh3/bwa))
- python packages ([pysam](http://pysam.readthedocs.io/en/latest/api.html), [pybedtools](https://daler.github.io/pybedtools/))
- R packages ([tsne](https://cran.r-project.org/web/packages/tsne/index.html), [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html), [densityClust](https://cran.r-project.org/web/packages/densityClust/index.html))

## Quick Install     
```bash
> git clone --depth=1 https://github.com/r3fang/snATAC.git
> cd snATAC
> chmod u+x bin/*
> export PATH=$PATH:./bin/
> ./bin/snATAC

Program:  snATAC (snATAC-seq analysis pipeline)
Version:  2.0.0
Contact:  Rongxin Fang <r3fang@ucsd.edu>
          Sebastian Preissl <spreissl@ucsd.edu>
          Bing Ren <biren@ucsd.edu>
          
Usage:    snATAC <command> [options]

Command:  pre           preprocessing
          filter        filter reads with unselected barcodes
          bmat          binary accessible matrix
          jacard        jaccard index matrix

Optional (under development): 
          decomplex     decomplex the fastq file
          bstat         simple statistics for barcode

Note: To use snATAC pipeline, you need to first decomplex barcode
combination and integrate barcodes to the beginning of the
read name in both R1 and R2 fastq files.
```

## Pipeline
**snATAC** has following steps:

1. Maping using [bwa](https://github.com/lh3/bwa) or [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (output.bam);
2. Pre-processing (output.bed/output.qc);
	* Alignment filteration
	* Speration of single cell
	* PCR duplicates removal
	* Mitochondrial reads removal
	* Adjusting position of Tn5 insertion
	* Create a master .bed file
	* Create a .qc file
3. Identify enrichment regions (output.txt);
	* Call peaks using MACS2 from ensemble data
4. Calculate barcode statistics;
	* Reads per barcode
	* Consecutive promoter coverage
	* Reads in peak ratio
4. Cell (barcode) selection (output.xgi);
	* Reads per barcode >= 1000
	* Consecutive promoter coverage > 5%
	* Reads in peak ratio >= 20%
	* **NOTE: above cutoff can vary singificant between different samples**
	* Filter potential doublets
5. Feature selection (output.ygi);
	* Filter top 5% peaks
	* Filter promoters
6. Binary accesibility matrix (output.mat);
7. Jaccard index matrix (output.jacard);
8. Dimentionality reduction (output.tsne);
9. Density-based clustering (output.cluster);

## A complete example

Step 0. download the sample data.     

```bash
> wget http://renlab.sdsc.edu/r3fang/snATAC/p56.rep1.R1.decomplex.fastq.gz
> wget http://renlab.sdsc.edu/r3fang/snATAC/p56.rep1.R2.decomplex.fastq.gz
> wget http://renlab.sdsc.edu/r3fang/snATAC/mm10.tar.gz
> tar -xvzf mm10.tar.gz
```

Step 1. Mapping using bwa (or bowtie2) in pair-end mode without any filteration.

```bash
> bwa mem -t 1 mm10.fa \
          p56.rep1.R1.decomplex.fastq.gz \
          p56.rep1.R2.decomplex.fastq.gz \
          | samtools view -bS - > p56.rep1.bam
```

Step 2. Pre-processing.

```bash
> ./bin/snATAC pre -t 5 -m 30 -f 2000 -e 75 \
                   -i p56.rep1.bam \
                   -o p56.rep1.bed.gz 2> p56.pre.log
```
This will output two files p56.pre.bed.gz and p56.pre.log contains basic quality control. Below is one example of p56.pre.log.

| p56.pre.log                         |          |
|:------------------------------------|:----------------|
| number of totol reads               |   210545337     |
| number of uniquely mapped reads     |   193513006     |
| number of properly paired reads     |   191507619     |
| number of chrM reads                |   2571611       |
| number of usable reads              |   188976583     |
| number of distinct reads            |   48747463      |
| estimated duplicate rate            |   0.74          |


Step 3. Identify peaks from ensemble signal

```bash
> macs2 callpeak -t p56.rep1.bed.gz \
                 -f BED -n p56.rep1 \
                 -g mm -p 0.05 \
                 --nomodel --shift 150 \
                 --keep-dup all    
                   
# extend peaks to 1kb from the summit
> bedtools flank -b 500 -g mm10.genome.size \
                 -i p56.rep1_summits.bed \
                 | sort -k1,1 -k2,2n - \
                 | bedtools merge -i - > p56.rep1.txt
```

Step 4. Cell (barcode) selection (output.xgi)

```bash
# count number of reads per barcode
> zcat p56.rep1.bed.gz | awk '{print $4}' \
               | sort \
               | uniq -c \
               | awk '{print $2, $1}' \
               | sort -k1,1 > p56.rep1.reads_per_cell
# consecutive promoter coverage 
> intersectBed -wa -wb \
               -a p56.rep1.bed.gz \
               -b mm10_consecutive_promoters.bed \
               | awk '{print $4, $8}' \
               | sort \
               | uniq \ 
               | awk '{print $1}' \
               | uniq -c \
               | awk '{print $2, $1}' \
               | sort -k1,1 > p56.rep1.promoter_cov 
# reads in peak ratio
> intersectBed -a p56.rep1.bed.gz -b p56.rep1.txt -u \
               | awk '{print $4}' \
               | sort \
               | uniq -c \
               | awk '{print $2, $1}' \
               | sort -k1,1 - > p56.rep1.reads_in_peak

```

Step 5. Cell selection (output.xgi)

```R
> R
##################################################################
# Reads per barcode >= 1000
# Consecutive promoter coverage > 3%
# Reads in peak ratio >= 20%
# NOTE: The cutoff can vary singificantly between different samples
##################################################################

consecutive_promoters <- read.table("mm10/mm10_consecutive_promoters.bed")
reads = read.table("p56.rep1.reads_per_cell")
promoters = read.table("p56.rep1.promoter_cov")
ratios = read.table("p56.rep1.reads_in_peak")
qc = reads; colnames(qc) <- c("barcode", "reads")
qc$promoter = 0; 
qc$read_in_peak = 0;
qc$promoter[match(promoters$V1, qc$barcode)] = 
	promoters$V2/nrow(consecutive_promoters)
qc$read_in_peak[match(ratios$V1, qc$barcode)] = ratios$V2
qc$ratio = qc$read_in_peak/qc$reads
idx <- which(qc$promoter > 0.03 & 
             qc$read_in_peak > 0.2 & 
             qc$reads > 1000)
qc_sel <- qc[idx,]

# among these cells, further filter barcodes with 
pvalues <- sapply(qc_sel$reads, function(x) 
           poisson.test(x, mean(qc_sel$reads), 
           alternative="greater")$p.value)
fdrs <- p.adjust(pvalues, "BH")
idx <- which(fdrs < 1e-2)

write.table(qc_sel[idx,1], file = "p56.rep1.xgi", append = FALSE, 
			  quote = FALSE, sep = "\t", eol = "\n", 
			  na = "NA", dec = ".", row.names = FALSE,
			  col.names = FALSE, qmethod = c("escape", "double"),
			  fileEncoding = "")
```

Step 5. Feature selection (output.xgi);

```R
> R
library(GenomicRanges)
peaks.df <- read.table("p56.rep1_peaks.narrowPeak")
# remove top 5% peaks
cutoff <- quantile((peaks.df$V5), probs = 0.95)
peaks.df <- peaks.df[which(peaks.df$V5 < cutoff),]
proms.df <- read.table("mm10/mm10.refSeq_promoter.bed")
peaks.gr <- GRanges(peaks.df[,1], IRanges(peaks.df[,2], peaks.df[,3]))
proms.gr <- GRanges(proms.df[,1], IRanges(proms.df[,2], proms.df[,3]))

peaks.sel.gr <- peaks.gr[-queryHits(findOverlaps(peaks.gr, proms.gr))]
peaks.sel.ex.gr <- resize(reduce(resize(peaks.sel.gr, 1000, 
                          fix="center")), 1000, fix="center")

peaks.sel.ex.df <- as.data.frame(peaks.sel.ex.gr)[,1:3]
write.table(peaks.sel.ex.df, file = "p56.rep1.ygi", 
			  append = FALSE, quote = FALSE, sep = "\t", 
			  eol = "\n", na = "NA", dec = ".", 
			  row.names = FALSE, col.names = FALSE, 
			  qmethod = c("escape", "double"),
			  fileEncoding = "")
```

Step 6. Generate binary accessibility matrix (**NOTE: this may require large RAM**)

```bash
> snATAC bmat -i p56.rep1.bed.gz \
              -x p56.rep1.xgi \
              -y p56.rep1.ygi \
              -o p56.rep1.mat
```

Step 7. Calculate jaccard index

```bash
# mannually count number of rows (1,464)
> wc -l p56.rep1.xgi
# mannually count number of columns (184,519)
> wc -l p56.rep1.ygi
# calculate jaccard index matrix
> snATAC jacard -i p56.rep1.mat -x 1464 -y 184519 -o p56.rep1.jacard
```

Step 8. Dimentionality reduction (R)
*(NOTE: tSNE is a heristic method, it is expected that different run can have slightly different result)*

```{R}
> R
library(tsne)
library(parallel)

data <- as.matrix(read.table("p56.rep1.jacard"))
diag(data) <- 0
b = tsne(data/sum(data), initial_config = NULL, 
		  k = 2, perplexity = 30, max_iter = 500, 
		  min_cost = 0, epoch_callback = NULL, 
		  whiten = TRUE, epoch=100)

plot(b, cex=0.7, xlab="tsne1", ylab="tsne2")

write.table(b, file = "p56.rep1.tsne", append = FALSE, 
            quote = FALSE, sep = "\t", eol = "\n", 
            na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
```
![plot](http://renlab.sdsc.edu/r3fang/snATAC/Rplots_tsne.pdf)

Step 9. Density-based clustering

```{R}
> R
library(densityClust)
MaxStep <- function(D){
	D_hat <- D
	n = nrow(D)
	for(k in 1:n){
		for(i in 1:(n-1)){
			for(j in (i+1):n){
				D_hat[i,j] <- min(D_hat[i,j], 
				max(D_hat[k,j], D_hat[i,k]))
			}
		}
	}
	return(D_hat)
}
find_center <- function(p, clust){
	centers <- data.frame()
	cols <- c()
	for(i in as.numeric(names(table(clust)))){
		ind <- which(clust== i)
		if(length(ind) > 10){
			centers <- rbind(centers, colMeans(p[ind,1:3]))
			cols <- c(cols, i)
		}
	}
	res <- cbind(centers, cols)
	colnames(res) <- c("x", "y", "z", "col")
	return(res)
}
DB_index <- function(D, CL){
	cl_uniq <- unique(CL)
	n <- length(cl_uniq)
	d_k <- rep(0, n)
	d_ij <- matrix(0, n, n)
	for(i in 1:n){
		ii <- which(CL == i);
		d_k[i] <- max(MaxStep(D[ii,ii]))
	}
	for(i in 1:(n-1)){
		for(j in (i+1):n){
			ii <- which(CL == i | CL == j)
			d_ij <- max(MaxStep(D[ii,ii]))
		}
	}
	return(min(d_ij)/max(d_k))
}

# perform cluster
points <- read.table("p56.rep1.tsne")

set.seed(10)
dis <- dist(points)
irisClust <- densityClust(dis, gaussian=TRUE)
# plot decision graph
plot(irisClust)
```
![Decision Graph](http://renlab.sdsc.edu/r3fang/snATAC/Rplots_decision_graph.pdf)

```R
# continue above code
rho_cutoff <- irisClust$dc
delta_cutoff <- 20
irisClust <- findClusters(irisClust, 
                          rho= rho_cutoff, 
                          delta= delta_cutoff)
clusters <- irisClust$cluster

plot(points, col=clusters, cex=0.5, pch=19)
dev.off()

write.table(data.frame(cluster), file = "p56.rep1.cluster",
			append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
```

![cluster](http://renlab.sdsc.edu/r3fang/snATAC/Rplots_cluster.pdf)

## Cite us
Preissl S.\*, Fang R.\*, Zhao Y., Raviram R., Zhang Y., Brandon C.S., Huang H., Gorkin D.U., Afzal V., Dickel D.E., Kuan S., Visel A., Pennacchio L.A., Zhang K., Ren B. **Single nucleus analysis of the chromatin landscape in mouse forebrain development**. bioRxiv 159137; doi: https://doi.org/10.1101/159137. (* contributed equally)




