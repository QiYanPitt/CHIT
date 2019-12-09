### generate haplotype_read_counts.$INDIVIDUAL.txt ###
ALL_SAMPLES_FILE=./geno_sample.lst
RNA_SAMPLES_FILE1=$1
for INDIVIDUAL in $(cat $RNA_SAMPLES_FILE1)
do
    python ./WASP/CHT/extract_haplotype_read_counts.py \
       --chrom ./chromInfo.hg19.txt \
       --snp_index results/snp_index.h5 \
       --snp_tab results/snp_tab.h5 \
       --haplotype results/haps.h5 \
       --geno_prob results/geno_probs.h5 \
       --samples $ALL_SAMPLES_FILE \
       --individual $INDIVIDUAL \
       --ref_as_counts results/ref_as_counts.$INDIVIDUAL.h5 \
       --alt_as_counts results/alt_as_counts.$INDIVIDUAL.h5 \
       --other_as_counts results/other_as_counts.$INDIVIDUAL.h5 \
       --read_counts results/read_counts.$INDIVIDUAL.h5 \
       ./gene_location/genes.gtf.wasp.input \
       > ./RNAseq/haplotypes/haplotype_read_counts.$INDIVIDUAL.txt
done

