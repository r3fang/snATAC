## Introduction
snATAC is a Ren-lab in-house pipeline for single-nucleus ATAC-seq [1,2] and snTHS-seq [3] analysis.

## Requirements
- [bwa](https://github.com/lh3/bwa) (recommend) or [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [samtools 1.2](http://www.htslib.org/doc/samtools.html)
- [Python 2.7](https://www.python.org/download/releases/2.7/)
- [bedtools](http://bedtools.readthedocs.io/en/latest/)
- [MACS2](https://github.com/taoliu/MACS) 
- [pysam](http://pysam.readthedocs.io/en/latest/api.html)
- [pybedtools](https://daler.github.io/pybedtools/)

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

Usage:    snATAC <command> [options]

Command:  pre           preprocessing
          filter        filter reads with unselected barcodes
          bmat          binary accessible matrix
          jacard        calculate jacard index matrix

Optional: decomplex     decomplex the fastq file
          bstat         simple statistics for barcode

Note: To use snATAC pipeline, you need to first decomplex barcode combination and integrate barcode information to the beginning of the read name in both R1 and R2 files. See FAQs for details.
```

## Pipeline
**snATAC** is made of following steps:

1. Maping using [bwa](https://github.com/lh3/bwa) or [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (output.bam);
2. Pre-processing (output.bed.gz);
	* Alignment filteration
	* Speration of single cell
	* PCR duplicates removal
	* Mitochondrial reads removal
	* Adjusting position of Tn5 insertion
	* Create a master file
	* Generate a .qc file
3. Identify peaks from ensemble signal using MACS2 (output.txt);
4. Cell (barcode) selection (output.xgi);
	* Number of reads
	* consecutive promoter coverage
	* Reads in peak ratio
5. Feature selection (output.ygi);
6. Binary accesibility matrix (output.mat);
7. Jaccard index matrix (output.jacard);
8. Dimentionality reduction (output.tsne);
9. Density-based clustering (output.cluster);
10. Identification of cluster-specific elements;
 
## A complete example

Step 0. download the sample data;     

```bash
> wget http://renlab.sdsc.edu/r3fang/snATAC/p56_R1.decomplex.fastq.gz
> wget http://renlab.sdsc.edu/r3fang/snATAC/p56_R2.decomplex.fastq.gz
```
Step 1. For snATAC-seq, mapping using bwa (or bowtie2) in pair-end mode without any filteration;

```bash
> bwa mem -t $NUM_CORES $MM10_BWA_INDEX \
          p56.rep1.R1.decomplex.fastq.gz \
          p56.rep1.R2.decomplex.fastq.gz \
          | samtools view -bS - > p56.rep1.bam
```
For snTHS-seq, mapping using bwa (or bowtie2) in single-end mode without any filteration;

```bash
> bwa mem -t $NUM_CORES $MM10_BWA_INDEX \
          p56.rep1.R1.decomplex.fastq.gz \
          | samtools view -bS - > p56.rep1.bam
```

Step 2. Pre-processing;

```bash
> ./bin/snATAC pre -t 5 -m 30 -f 2000 -e 75 \
                   -i p56.rep1.bam \
                   -o p56.rep1.pre.bed.gz 2> p56.pre.log
```
This will output two files p56.pre.bed.gz and p56.pre.log. p56.pre.bed.gz.qc contains basic quality control. Below is one example of p56.pre.bed.gz.qc

| Item                               | Number          |
|:------------------------------------|:----------------|
| number of totol reads               |   210545337     |
| number of uniquely mapped reads     |   193513006     |
| number of properly paired reads     |   191507619     |
| number of chrM reads                |   2571611       |
| number of usable reads              |   188976583     |
| number of distinct reads            |   48747463      |
| estimated duplicate rate            |   0.742044      |


Step 3. Identify peaks from ensemble signal

```bash
> macs2 callpeak   -t p56.rep1.bed.gz \
                   -f BED -n p56.rep1 \
                   -g mm -p 0.05 \
                   --nomodel --shift 150 \
                   --keep-dup all    
                   
# extend peaks to 1kb from the summit
> bedtools flank   -b 500 -g mm10.genome \
                   -i p56.rep1_summits.bed \
                   | sort -k1,1 -k2,2n - \
                   | bedtools merge -i - > p56.rep1.txt
```

Step 4. Cell (barcode) selection (.xgi);

```bash
# count number of reads per cell
> zcat p56.rep1.pre.bed.gz \
        | awk '{print $3}' \
        | sort | uniq -c > p56.rep1.number_of_reads

# consecutive promoter coverage 
> 
# reads in peak ratio
> 
# select cells
> 
```

Step 5. Feature selection (.ygi);

```bash
# calculate feature coverage
> 
# filter promoter regions and high-coverage peaks
> intersectBed -v  -a p56.features.bed \
                   -b refSeq_promoter.bed > p56.features.distal.bed 
```

Step 6. Generate binary accessibility matrix

```bash
> snATAC bmat -i p56.rep1.bed.gz \
              -x p56.rep1.xgi \
              -y p56.rep1.ygi \
              -o p56.rep1.mat
```

Step 7. Calculate jaccard index

```bash
# mannually count number of rows
> wc -l p56.rep1.xgi
# mannually count number of columns
> wc -l p56.rep1.ygi
# calculate jaccard index matrix
> snATAC_jacard -i p56.rep1.mat -x 1500 -y 7856 -o p56.rep1.jacard
```

Step 8. Dimentionality reduction

```{R}
> R
library(tsne)
library(parallel)

data <- as.matrix(read.table("p56.rep1.jacard"))
diag(data) <- 0
b = tsne(data/sum(data), initial_config = NULL, k = 2, perplexity = 30, max_iter = 2000, min_cost = 0, epoch_callback = NULL, whiten = TRUE, epoch=100)

write.table(b, file = "p56.rep1.tsne", append = FALSE, 
            quote = FALSE, sep = "\t", eol = "\n", 
            na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
```

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
				D_hat[i,j] <- min(D_hat[i,j], max(D_hat[k,j], D_hat[i,k]))
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
irisClust <- findClusters(irisClust, rho=20, delta=50)
cluster <- irisClust$cluster

write.table(data.frame(cluster), file = args[3],
			append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
```
