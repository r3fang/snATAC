##Get Started     
```bash
$ git clone https://github.com/r3fang/scATAC.git
$ cd scATAC
$ chmod u+x bin/*
$ export PATH=$PATH:./bin/ #(add this path to .bash_profile permanently)
$ scATAC_debarcode -a data/demo_ATAC.I1.gz -b data/demo_ATAC.I2.gz -c data/demo_ATAC.R1.gz \
	     | gzip - > demo_ATAC.decomplex.R1.fastq.gz
$ scATAC_debarcode -a data/demo_ATAC.I1.gz -b data/demo_ATAC.I2.gz -c data/demo_ATAC.R2.gz \
	     | gzip - > demo_ATAC.decomplex.R2.fastq.gz
$ scATAC -t 2 \
         -f demo_ATAC.decomplex.R1.fastq.gz \
         -r demo_ATAC.decomplex.R2.fastq.gz \
		 -b ./barcodes \
		 -d 2 \
		 -p Picard/MarkDuplicates.jar \
		 -n demo_atac \
		 -g hg19.fa \
		 -m 500
```
##Depedency
- [bwa](https://github.com/lh3/bwa)
- [samtools 1.2+](http://www.htslib.org/doc/samtools.html)
- [Python 2.7+](https://www.python.org/download/releases/2.7/)

##Introduction

**scATAC** is an in-house Bioinformatics pipeline for analyzing multiplex single-cell ATAC-seq data.

```
$ scATAC

Program: scATAC (Multiplex single-cell ATAC-seq analysis pipeline by Ren Lab)
Version: 06.10.2016
Contact: Rongxin Fang <r3fang@ucsd.edu>
Ren Lab: http://bioinformatics-renlab.ucsd.edu/rentrac/

usage: scATAC [-h] [-t THREADS] [-f R1] [-r R2] [-b BARCODE_DIR] [-d MAX_BARCODE_MISMATCH] [-p MarkDuplicates.jar] [-n PREFIX] [-g BWA_GENOME] [-m MIN_READ]

Example:
scATAC -t 2 -f demo_R1.fastq.gz -r demo_R2.fastq.gz -b ./barcodes -d 2 -p Picard/MarkDuplicates.jar -n demo -g hg19.fa -m 500

Options:
	-h, --help                  show this help message and exit.
	-t  THREADS                 threads [1].
	-f  R1                      fastq.gz file that contains forward reads (only fastq.gz allowed).
	-r  R2                      fastq.gz file that contains reverse reads (only fastq.gz allowed).
	-b  BARCODE_DIR             folder that contains r7_ATAC, i7_ATAC, i5_ATAC and r5_ATAC barcode.
	-d  MAX_BARCODE_MISMATCH    max barcode mismatch allowed for barcode error correction [2].
	-p  MARK_DUPLICATE          path to picard MarkDuplicates.jar [MarkDuplicates.jar].
	-n  PREFIX                  prefix of output files.
	-g  BWA_GENOME              BWA indexed reference genome.
	-m  MIN_READ                cells with reads less than MIN_READ will be filtered [500].

Note: To use scATAC, you need to first decomplex barcode combination and integrate barcode
      information as the beginning of the read name in R1 and R2 files.
      This can be done by command 'scATAC_debarcode':
      'scATAC_debarcode -a I1.fastq.gz -b I2.fastq.gz -c R1.fastq.gz | gzip - > R1.decomplex.fastq.gz'
      'scATAC_debarcode -a I1.fastq.gz -b I2.fastq.gz -c R2.fastq.gz | gzip - > R2.decomplex.fastq.gz'
```

##Pipeline

**scATAC** is made of following steps:

0. decomplex scATAC-seq data by scATAC_debarcode [OPTIONAL];
1. map using bwa followed by filtering reads with MAPQ < 10;
2. correct barcode error caused by sequencing error by allowing certain number of mismatches [2];
3. split reads to individual cells based on the barcode combination;
4. remove PCR duplication for each cell;
6. merge reads from different cells;
7. generate barcode frequency table;
8. filter cells with reads counts less than given number [500];
9. summerize and generate a log file;

##Output
**scATAC** generates two files '.log' and '.bam'. 
'.bam' is the final file that contains all usable reads and '.log' includes data metrics.

##FAQ


 1. **How to remove mitochondrial reads from BAM files?**

 ```bash
 $samtools index out.bam
 $samtools idxstats out.bam | cut -f 1 | grep -v chrM | xargs samtools view -b - > out.filtered.bam
 ```

 2. **How to get barcode fequency if reads mapped to concatenated genome?**
 
 ```bash
 # generate barcode frequency for mm9
 samtools view out.bam | awk '$3 ~ /mm9/ {split($1,a,":"); print a[1]}' \
 	| sort | uniq -c | awk '{print $2, $1}' | sort -k2rn - > out.mm9.barcode_freq.txt
 # generate barcode frequency for hg19
 samtools view out.bam | awk '$3 ~ /hg19/ {split($1,a,":"); print a[1]}' \
 	| sort | uniq -c | awk '{print $2, $1}' | sort -k2rn - > out.hg19.barcode_freq.txt
 ```




