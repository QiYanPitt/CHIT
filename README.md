# Allele-specific method for testing the association between gene expression and phenotype-genotype interaction

## Introduction
The combined haplotype interaction test (CHIT) tests the association between gene expressions and phenotype-genotype interactions by modeling the total read counts and allele-specific reads in a target gene region. Each test is performed based on a test SNP-feature region pair that is defined by users. Thus, this method can test both cis- and trans-regulatory effect.

The detailed pipeline is described in the README file, including:
* Preprocess phased SNP files;
* Generate SNP-gene pairs file;
* Use both SNP and bam (e.g., RNA-seq) files to remap bam, done by WASP (WASP: allele-specific software for robust molecular quantitative trait locus discovery [PMID: 26366987]);
* Count total read counts and allele-specific reads in a target gene region, done by WASP;
* Estimate overdispersion parameters, done by WASP;
* Run CHIT;
* Get concordant results from different initial values.

## Dependencies
The dependencies are the same as WASP that is described in ./WASP/README.md

## Usage of 
for i in {1..22}; do
    CHIT_IN_FILE=./RNAseq/cht_input_file_atopy_ext.ranphe.txt.$i
    nohup python ./CHIT/CHIT.py \
        --bnb_disp ./RNAseq/cht_bnb_coef_atopy.txt \
        --as_disp ./RNAseq/cht_as_coef_atopy.txt \
        --cov_file ./pca/cov.matrix.txt \
        $CHIT_IN_FILE ./RNAseq/output/chit_atopy_results.txt.$i
done

## Options
