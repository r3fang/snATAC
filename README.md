##Get Started     
```
$ git clone https://github.com/r3fang/scATAC.git
$ cd scATAC
$ chmod u+x bin/
$ export PATH=$PATH:bin/
$ scATAC_debarcode -a data/demo_ATAC.I1.gz -b data/demo_ATAC.I2.gz -c data/demo_ATAC.R1.gz | gzip - > demo_ATAC.decomplex.R1.fastq.gz
$ scATAC_debarcode -a data/demo_ATAC.I1.gz -b data/demo_ATAC.I2.gz -c data/demo_ATAC.R2.gz | gzip - > demo_ATAC.decomplex.R2.fastq.gz
$ scATAC -t 2 -f demo_ATAC.decomplex.R1.fastq.gz -r demo_ATAC.decomplex.R2.fastq.gz -b ./barcodes -d 2 -p Picard/MarkDuplicates.jar -n demo_ATAC -g hg19.fa -m 500 
```

##Introduction

**scATAC** is a high-efficient in-house Bioinformatics pipeline designed for analyzing dual-barcode signle-cell ATAC-seq data.

```
$ scATAC

Program: scATAC (dual-barcode single cell ATAC-seq analysis pipeline by Ren Lab)
Version: 06.10.2016
Contact: Rongxin Fang <r3fang@ucsd.edu>
Ren Lab: http://bioinformatics-renlab.ucsd.edu/rentrac/

Step 0. decomplex scATAC-seq data [OPTIONA];
Step 1. map and filter poorly mapped reads;
Step 2. correct barcode error by allowing given number of mismatches;
Step 3. split reads to individual cells based on the barcode combination;
Step 4. remove PCR duplication for each cell;
Step 6. merge reads from different cells;
Step 7. generate barcode frequency table;
Step 8. filter cells with reads counts less than given cutoff;
Step 9. summerize and generate a log file;

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
	-m  MIN_READ                cells with reads less than MIN_READ will be filtered.

Note: To use scATAC, you need to first decomplex barcode combination and integrate barcode
      information as the beginning of the read name in R1 and R2 files.
      This can be done by command 'scATAC_debarcode':
      'scATAC_debarcode -a I1.fastq.gz -b I2.fastq.gz -c R1.fastq.gz | gzip - > R1.decomplex.fastq.gz'
      'scATAC_debarcode -a I1.fastq.gz -b I2.fastq.gz -c R2.fastq.gz | gzip - > R2.decomplex.fastq.gz'
```
