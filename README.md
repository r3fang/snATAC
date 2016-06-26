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
##Requirements

- [bwa](https://github.com/lh3/bwa)
- [samtools 1.2+](http://www.htslib.org/doc/samtools.html)
- [Python 2.7+](https://www.python.org/download/releases/2.7/)
- [bedtools](http://bedtools.readthedocs.io/en/latest/)
- [MACS2](https://github.com/taoliu/MACS) (OPTIONAL: needed only for peak calling)
- [Numpy](http://www.numpy.org/) [OPTIONAL: needed only for generating accessible binary matrix]
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
 
##Complete Example

```bash
 $ git clone https://github.com/r3fang/scATAC.git
 $ cd scATAC
 $ chmod u+x bin/*
 $ export PATH=$PATH:./bin/ #(add the absolute path to .bash_profile permanently)

 # download decomplexed sample data (human GM12878 and mouse ES mixture) from Cusanovich, Science, 2015
 $ wget http://enhancer.sdsc.edu/r3fang/Cusanovich_2015/SRR1947691_1.fastq.gz
 $ wget http://enhancer.sdsc.edu/r3fang/Cusanovich_2015/SRR1947691_2.fastq.gz
 
 # download BWA-indexed reference genome hg19-mm9 concatenated genome
 $ wget http://enhancer.sdsc.edu/r3fang/Cusanovich_2015/hg19_mm9_genome.tar.gz
 $ tar -xvzf hg19_mm9_genome.tar.gz
 
 # map and analyze
 $ scATAC -t 5 -f SRR1947691_1.fastq.gz -r SRR1947691_2.fastq.gz -b barcodes/ -d 1 \
 -p Picard/MarkDuplicates.jar -n SRR1947691 -g hg19_mm9_genome/genome.fa -m 500 > SRR1947691.log
 
 # remove mitochondrial reads
 $ samtools index SRR1947691.bam
 $ samtools idxstats SRR1947691.bam | cut -f 1 | grep -v chrM \
   | xargs samtools view -bS SRR1947691.bam > SRR1947691.filtered.bam
 
 # split reads into human and mouse
 $ samtools view -h SRR1947691.filtered.bam \
   | scATAC_split_genome hg19 \
   | samtools view -bS - >  SRR1947691.filtered.hg19.bam

 $ samtools view -h SRR1947691.filtered.bam \
   | scATAC_split_genome mm9 \
   | samtools view -bS - >  SRR1947691.filtered.mm9.bam

 # convert bulk bam to bigWig file
 $ fetchChromSizes hg19 > hg19.chrom.sizes
 $ fetchChromSizes mm9 > mm9.chrom.sizes
 
 $ bamToBed -i SRR1947691.filtered.hg19.bam \
   | slopBed -s -l 0 -r 300 -i stdin -g hg19.chrom.sizes \
   | bedtools genomecov -g hg19.chrom.sizes -i stdin -bg \
   | sort -k1,1 -k2,2n - | wigToBigWig stdin hg19.chrom.sizes SRR1947691.filtered.hg19.bw 
 
 $ bamToBed -i SRR1947691.filtered.mm9.bam \
   | slopBed -s -l 0 -r 300 -i stdin -g mm9.chrom.sizes \
   | bedtools genomecov -g mm9.chrom.sizes -i stdin -bg \
   | sort -k1,1 -k2,2n - | wigToBigWig stdin mm9.chrom.sizes SRR1947691.filtered.mm9.bw 
 
 # generate barcode frequency
 $ samtools view SRR1947691.filtered.hg19.bam | awk '{split($1,a,":"); print a[1]}' \
   | sort | uniq -c | awk '{print $2, $1}' | sort -k2rn - > SRR1947691.filtered.hg19.barcode_freq.txt 
 
 $ samtools view SRR1947691.filtered.mm9.bam | awk '{split($1,a,":"); print a[1]}' \
   | sort | uniq -c | awk '{print $2, $1}' | sort -k2rn - > SRR1947691.filtered.mm9.barcode_freq.txt
 
 # peak calling using MACS2
 $ macs2 callpeak -t SRR1947691.filtered.hg19.bam -f BAM -g hs \
   --outdir SRR1947691.filtered.hg19 -n SRR1947691.filtered.hg19 --extsize 300 -q 0.1
 $ macs2 callpeak -t SRR1947691.filtered.mm9.bam -f BAM -g mm \
   --outdir SRR1947691.filtered.mm9 -n SRR1947691.filtered.mm9 --extsize 300 -q 0.1   
```

## FAQ

 1. **How to remove reads from .bam file whose barcodes are not selected?**     
 First, you need to provide a barcodes.sel.txt file that has barcodes you want to keep seperated by space or lines. Then run following command

 ```bash
 # filter reads whose barcodes are not selected
 $ samtools view -h input.bam \
   | scATAC_rm_reads_by_barcodes barcodes.sel.txt - | samtools view -bS - >  out.filtered.bam 
 ```

 2. **How to generate accessible binary matrix?**     
 First, you need to provide (1) a .bed file peaks.bed has all the inquiry regions and (2) a barcode.sel.txt file has a list of inquiry barcodes. Next, be sure that .bed file is a valid `awk '{if($3 <= $2) print }' peaks.bed`, it is valid if nothing printed on the screen, otherwise change manuanlly or `awk '{if($3 <= $2) printf "%s\t%d\t%d\t%s\n", $1, $3, $2, $4; else printf "%s\t%d\t%d\t%s\n", $1, $2, $3, $4 }' peaks.bed > peaks.valid.bed`. Finally you can generate .mat, .xgi, .ygi by
 '''bash
 $ awk '{printf "%s\t%d\t%d\t%s\n", $1, $2, $3, $4}' peaks.bed \
   | intersectBed -wa -wb -a stdin -b input.bam \
   | awk '{print $1, $2, $3, $4, $8}' \
   | scATAC_get_binary_mat peaks.bed barcode.sel.txt prefix 
 '''



