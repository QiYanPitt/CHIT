# Allele-specific method for testing the association between gene expression and phenotype-genotype interaction

## Introduction
The combined haplotype interaction test (CHIT) tests the association between gene expressions and phenotype-genotype interactions by modeling the total read counts and allele-specific reads in a target gene region. Each test is performed based on a test SNP-feature region pair that is defined by users. Thus, this method can test both cis- and trans-regulatory effect.

The detailed pipeline is described in the [README](./README) file, including:
* Preprocess phased SNP files;
* Generate SNP-gene pairs file;
* Use both SNP and bam (e.g., RNA-seq) files to remap bam, done by WASP (WASP: allele-specific software for robust molecular quantitative trait locus discovery [PMID: 26366987]);
* Count total read counts and allele-specific reads in a target gene region, done by WASP;
* Estimate overdispersion parameters, done by WASP;
* Run CHIT;
* Get concordant results from different initial values.

## Dependencies
The dependencies are the same as WASP that is described in [./WASP/README.md](./WASP/README.md)

## CHIT usage
    for i in {1..22}; do
        CHIT_IN_FILE=./RNAseq/cht_input_file_atopy_ext.ranphe.txt.$i
            python ./CHIT/CHIT.py \
            --bnb_disp ./RNAseq/cht_bnb_coef_atopy.txt \
            --as_disp ./RNAseq/cht_as_coef_atopy.txt \
            --cov_file ./pca/cov.matrix.txt \
            $CHIT_IN_FILE ./RNAseq/output/chit_atopy_results.txt.$i
    done

## CHIT options
    * "-a", "--as_only": only perform the allele-specific part (Beta Binomial) part of the test;
    * "-d", "--bnb_only": only perform the association (Beta Negative Binomial) part of the test;
    * "--cov_file": file containing covariates to include in the model;
    * "-b", "--bnb_disp": file containing depth (Beta Negative Binomial) dispersion parameters;
    * "-o", "--as_disp": file containing allele-specific (Beta Binomial) dispersion parameters;
    * "-s", "--shuffle": permute genotypes;
    * "-m", "--min_as_counts": only perform test when total number of allele-specific read counts across individuals > MIN_COUNTS;
    * "-v", "--verbose": print extra information;
    * "--benchmark": write information about time test is takes, number of optimization functions, etc.;
    * "-i", "--initial": initial value for parameters alpha and beta;
    * "-r", "--seed": seed for random shuffle;
    * infile_list;
    * out_file;
