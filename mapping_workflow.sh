#!/bin/bash
DATA_DIR=./
RNA_SAMPLES_FILE1=$1
ALL_SAMPLES_FILE=./geno_sample.lst
genome_STAR_index=/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/STARIndex #need to download from UCSC
threads=1

### Pull out reads that need to be remapped to check for bias ###
# Use the -p option for paired-end reads.
mkdir ttt
for INDIVIDUAL in $(cat $RNA_SAMPLES_FILE1)
do
  echo ${INDIVIDUAL} > ttt/${INDIVIDUAL}.lst
  RNA_SAMPLES_FILE=ttt/${INDIVIDUAL}.lst
  python ./WASP/mapping/find_intersecting_snps.py \
       --is_paired_end \
       --output_dir $DATA_DIR/RNAseq  \
       --snp_index $DATA_DIR/results/snp_index.h5 \
       --snp_tab $DATA_DIR/results/snp_tab.h5 \
       --haplotype $DATA_DIR/results/haps.h5 \
       --samples $RNA_SAMPLES_FILE \
       ./bam/${INDIVIDUAL}.bam &
done
wait
rm -r ttt

for INDIVIDUAL in $(cat $RNA_SAMPLES_FILE1); do
  rm $DATA_DIR/RNAseq/${INDIVIDUAL}.remap.single.fq.gz
  rm $DATA_DIR/RNAseq/${INDIVIDUAL}.sort.bam
done

### Remap the reads, using same the program and options as before ###
mkdir $DATA_DIR/tmp
STAR --genomeDir ${genome_STAR_index} --genomeLoad LoadAndExit --outFileNamePrefix $DATA_DIR/tmp/
rm -rdf $DATA_DIR/tmp
for INDIVIDUAL in $(cat $RNA_SAMPLES_FILE1)
do
  STAR --genomeDir ${genome_STAR_index} --outSAMunmapped Within  --outFilterType BySJout  --outSAMattributes NH HI AS NM MD  --outFilterMultimapNmax 20  --outFilterMismatchNmax 999  --outFilterMismatchNoverLmax 0.04  --alignIntronMin 20  --alignIntronMax 1000000  --alignMatesGapMax 1000000  --alignSJoverhangMin 8  --alignSJDBoverhangMin 1  --sjdbScore 1  --runThreadN ${threads}  --genomeLoad LoadAndKeep  --outSAMtype BAM Unsorted  --quantMode TranscriptomeSAM  --outSAMheaderHD \@HD VN:1.4 SO:unsorted  --outFileNamePrefix $DATA_DIR/RNAseq/${INDIVIDUAL}.remap  --readFilesCommand zcat  --readFilesIn $DATA_DIR/RNAseq/${INDIVIDUAL}.remap.fq*.gz &
done
wait
rm $DATA_DIR/RNAseq/*SJ.out.tab $DATA_DIR/RNAseq/*Log.*
echo "=== analysis finished, release genome from shared memory ==="
mkdir $DATA_DIR/tmp
STAR --genomeDir ${genome_STAR_index} --genomeLoad Remove --outFileNamePrefix $DATA_DIR/tmp/
rm -rdf $DATA_DIR/tmp

### Use filter_remapped_reads.py to create filtered list of reads that correctly ###
# remap to same position
for INDIVIDUAL in $(cat $RNA_SAMPLES_FILE1); do
  python ./WASP/mapping/filter_remapped_reads.py \
      $DATA_DIR/RNAseq/${INDIVIDUAL}.to.remap.bam \
      $DATA_DIR/RNAseq/${INDIVIDUAL}.remapAligned.out.bam \
      $DATA_DIR/RNAseq/${INDIVIDUAL}.remap.keep.bam &
  ### Create a merged BAM containing [1] reads that did
  ### not need remapping [2] filtered remapped reads
done
wait

for INDIVIDUAL in $(cat $RNA_SAMPLES_FILE1); do
  rm $DATA_DIR/RNAseq/${INDIVIDUAL}.to.remap.bam
  rm $DATA_DIR/RNAseq/${INDIVIDUAL}.remapAligned.toTranscriptome.out.bam
  rm $DATA_DIR/RNAseq/${INDIVIDUAL}.remapAligned.out.bam
  samtools merge $DATA_DIR/RNAseq/${INDIVIDUAL}.keep.merged.bam \
           $DATA_DIR/RNAseq/${INDIVIDUAL}.keep.bam $DATA_DIR/RNAseq/${INDIVIDUAL}.remap.keep.bam &
done
wait

for INDIVIDUAL in $(cat $RNA_SAMPLES_FILE1); do
  rm $DATA_DIR/RNAseq/${INDIVIDUAL}.keep.bam
  rm $DATA_DIR/RNAseq/${INDIVIDUAL}.remap.keep.bam
  samtools sort $DATA_DIR/RNAseq/${INDIVIDUAL}.keep.merged.bam -o $DATA_DIR/RNAseq/${INDIVIDUAL}.keep.merged.sorted.bam &
done
wait

for INDIVIDUAL in $(cat $RNA_SAMPLES_FILE1); do
  rm $DATA_DIR/RNAseq/${INDIVIDUAL}.keep.merged.bam
  samtools index $DATA_DIR/RNAseq/${INDIVIDUAL}.keep.merged.sorted.bam
  ### Filter out duplicate reads. Use rmdup_pe.py for paired-end reads, rmdup.py for single-end reads.
  python ./WASP/mapping/rmdup_pe.py $DATA_DIR/RNAseq/${INDIVIDUAL}.keep.merged.sorted.bam $DATA_DIR/RNAseq/${INDIVIDUAL}.keep.rmdup.merged.sorted.bam &
done
wait

for INDIVIDUAL in $(cat $RNA_SAMPLES_FILE1); do
  rm $DATA_DIR/RNAseq/${INDIVIDUAL}.keep.merged.sorted.bam
  rm $DATA_DIR/RNAseq/${INDIVIDUAL}.keep.merged.sorted.bam.bai
  samtools sort $DATA_DIR/RNAseq/${INDIVIDUAL}.keep.rmdup.merged.sorted.bam -o $DATA_DIR/RNAseq/${INDIVIDUAL}.keep.rmdup.merged.sorted2.bam &
done
wait

for INDIVIDUAL in $(cat $RNA_SAMPLES_FILE1); do
  rm $DATA_DIR/RNAseq/${INDIVIDUAL}.keep.rmdup.merged.sorted.bam
  samtools index $DATA_DIR/RNAseq/${INDIVIDUAL}.keep.rmdup.merged.sorted2.bam &
done 
wait

### generate ref_as_counts.$INDIVIDUAL.h5, alt_as_counts.$INDIVIDUAL.h5, other_as_counts.$INDIVIDUAL.h5 and read_counts.$INDIVIDUAL.h5 ###
for INDIVIDUAL in $(cat $RNA_SAMPLES_FILE1)
do
    echo $INDIVIDUAL
    python ./WASP/CHT/bam2h5.py --chrom ./chromInfo.hg19.txt \
              --snp_index results/snp_index.h5 \
              --snp_tab results/snp_tab.h5 \
              --haplotype results/haps.h5 \
              --samples $ALL_SAMPLES_FILE \
              --individual $INDIVIDUAL \
              --ref_as_counts results/ref_as_counts.$INDIVIDUAL.h5 \
              --alt_as_counts results/alt_as_counts.$INDIVIDUAL.h5 \
              --other_as_counts results/other_as_counts.$INDIVIDUAL.h5 \
              --read_counts results/read_counts.$INDIVIDUAL.h5 \
              $DATA_DIR/RNAseq/$INDIVIDUAL.keep.rmdup.merged.sorted2.bam &
done 
wait

for INDIVIDUAL in $(cat $RNA_SAMPLES_FILE1); do
    rm $DATA_DIR/RNAseq/$INDIVIDUAL.keep.rmdup.merged.sorted2.bam
    rm $DATA_DIR/RNAseq/$INDIVIDUAL.keep.rmdup.merged.sorted2.bam.bai
done

