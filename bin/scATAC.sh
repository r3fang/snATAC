#!/bin/bash

# PART I check weather softwares installed
command -v bwa >/dev/null 2>&1 || { echo >&2 "scATAC requires bwa but it is not installed.  Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "scATAC requires samtools but it is not installed.  Aborting."; exit 1; }



PREFIX='Undetermined_S0_L001' # make sure your imput fastq files is $PREFIX.R1.fastq.gz and $PREFIX.R2.fastq.gz
GENOME='/oasis/tscc/scratch/r3fang/data/mixture_hg19_mm9/Sequence/WholeGenomeFasta/genome.fa'
GENOME_SIZE='/oasis/tscc/scratch/r3fang/data/mixture_hg19_mm9/Sequence/WholeGenomeFasta/genome.fai'
TMP_FOLDER='Undetermined_S0_L001.tmp'
scATAC_barcode_err_correct='/oasis/tscc/scratch/r3fang/collaboration/Seb_04_26_2016/Mixtures_published_download/HL60_GM12878/scATAC_barcode_err_correct.py'
scATAC_decell='/oasis/tscc/scratch/r3fang/collaboration/Seb_04_26_2016/Mixtures_published_download/HL60_GM12878/scATAC_decell.py'
mark_duplcate='/oasis/tscc/scratch/r3fang/github/Hicmap_Ren/Picard/MarkDuplicates.jar'
scATAC_rm_cell_with_low_cov='/oasis/tscc/scratch/r3fang/collaboration/Seb_04_26_2016/Mixtures_published_download/HL60_GM12878/scATAC_rm_cell_with_low_cov.py'
# echo 'Step 1. map and filter reads with poor mappability'
# bwa mem -M -t 5 $GENOME $PREFIX.R1.fastq.gz $PREFIX.R2.fastq.gz \
# | samtools view -q 10 -bS - > $PREFIX.umap.bam

echo 'Step 2. correct barcode allowing 2 mismatches'
samtools view -h $PREFIX.umap.bam \
| python $scATAC_barcode_err_correct 2 \
| samtools view -bS - > $PREFIX.umap.2mm.bam

echo 'Step 3. sort reads based on name'
samtools sort -n -m 1G $PREFIX.umap.2mm.bam $PREFIX.umap.2mm.nsorted

echo 'Step 4. seperate reads into each cell' 
mkdir $TMP_FOLDER # create a tmp folder
samtools view $PREFIX.umap.2mm.nsorted.bam \
| python $scATAC_decell $TMP_FOLDER -

echo 'Step 5. remove PCR duplication for each cell'
samtools view -H $PREFIX.umap.2mm.nsorted.bam > $PREFIX.umap.2mm.nsorted.header
for line in `ls $TMP_FOLDER | grep .sam`
do
	barcode="${line%.*}"
	cat $PREFIX.umap.2mm.nsorted.header $TMP_FOLDER/$line | samtools view -bS - | samtools sort - $TMP_FOLDER/$barcode.sorted
	java -Xmx2G -jar $mark_duplcate INPUT=$TMP_FOLDER/$barcode.sorted.bam OUTPUT=$TMP_FOLDER/$barcode.sorted.filtered.bam ASSUME_SORTED=true \
		REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE=$TMP_FOLDER/metrics.$barcode.txt TMP_DIR=$TMP_FOLDER/$barcode\_tmp
	rm $TMP_FOLDER/$barcode.sorted.bam
	rm -r $TMP_FOLDER/$barcode\_tmp
done

echo 'Step 6. merge cells'
samtools cat -o $PREFIX.umap.2mm.nsorted.nodup.bam $TMP_FOLDER/*.sorted.filtered.bam

echo 'Step 7. remove intermedia files/folders'
rm $PREFIX.umap.2mm.nsorted.header
rm -r $TMP_FOLDER

echo 'Step 8. sort by genomic coordinates'
samtools sort -m 1G $PREFIX.umap.2mm.nsorted.nodup.bam $PREFIX.umap.2mm.nsorted.nodup.gsorted

echo 'Step 9. generate barcode frequency table'
samtools view $PREFIX.umap.2mm.bam \
| awk '{split($1,a,":"); print a[1]}' | sort | uniq -c | awk '{print $2, $1}' \
| sort -k2rn - > $PREFIX.umap.2mm.stat &

samtools view $PREFIX.umap.2mm.nsorted.nodup.gsorted.bam \
| awk '{split($1,a,":"); print a[1]}' | sort | uniq -c | awk '{print $2, $1}' \
| sort -k2rn - > $PREFIX.umap.2mm.nsorted.nodup.gsorted.stat 

echo 'Step 10. filter cells with low reads (less than 500 reads)'
samtools view -h $PREFIX.umap.2mm.nsorted.nodup.gsorted.bam \
| python $scATAC_rm_cell_with_low_cov 500 $PREFIX.umap.2mm.nsorted.nodup.gsorted.stat - \
| samtools view -bS - > $PREFIX.umap.2mm.nsorted.nodup.gsorted.filtered.bam

#echo 'Step 11. convert bam file to bigwig for visualization'
#bamToBed -i $PREFIX.umap.2mm.nsorted.nodup.gsorted.bam \
#| slopBed -s -l 0 -r 300 -i stdin -g hg19.genome \
#| genomeCoverageBed -g hg19.genome -i stdin -bg \
#| wigToBigWig stdin hg19.genome SRR1947693.umap.2mm.nsorted.uniq.filtered.gsorted.bw
