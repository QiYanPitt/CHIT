# Copyright 2013 Graham McVicker and Bryce van de Geijn
import sys
import os
import math
import time
import gzip
import argparse

from scipy.optimize import *
from scipy.special import gammaln
from scipy.special import betaln
import scipy.stats

import numpy as np

import random
import util 
import test_parse_test_snp
import parse_options
import rescale_totals
import ll
import write_results
import write_header
import read_as_sigmas
import read_bnb_sigmas

OPTIMIZER="Nelder-Mead"
#OPTIMIZER="BFGS"
options = parse_options.parse_options()
options.dup_snp_warn = True
cov_matrix = []

def write_empty_result(outfile, snpinfo):
    """Write all zeros in the even that the test failed"""
    outfile.write("\t".join([snpinfo[0][0], snpinfo[0][1], "0", "0",
                             "0", "NA", "0", "0", "0", "0",
                             "0", "0", "0"]) + '\n')

def open_input_files(in_filename):
    # read file that contains list of input files
    in_file = open(in_filename)
    infiles = []
    for line in in_file:
        # open each input file and read first line
        filename = line.rstrip().split(" ")[0]
        if util.is_gzipped(filename):
            f = gzip.open(filename, "rt")
        else:
            f = open(filename, "r")
        # skip header
        f.readline()
        infiles.append(f)
    in_file.close()
    return infiles

def extract_disease_status(in_filename):
    in_file = open(in_filename)
    disease_stat = []
    for line in in_file:
        stat = line.rstrip().split(" ")[1]
        disease_stat.append(stat)
    in_file.close()
    return disease_stat

disease_stat = extract_disease_status(options.infile_list)
infiles = open_input_files(options.infile_list)
outfile = open(options.out_file, 'w')
write_header.write_header(outfile)

snpinfo = []
for f in infiles:
    snpinfo.append(f.readline().strip().split())

row_count = 0
finished=False

# read dispersion parameters for each individual
bnb_sigmas = read_bnb_sigmas.read_bnb_sigmas(options, infiles)
as_sigmas = read_as_sigmas.read_as_sigmas(options, infiles)

###########################
if options.cov_file:
    # convert covariate matrix to pca matrix
    cov_matrix = ll.load_covariates(options.cov_file)
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    cov_matrix = StandardScaler().fit_transform(cov_matrix)
    pca = PCA(n_components=cov_matrix.shape[1])
    cov_matrix = pca.fit_transform(cov_matrix)
    num_covs = cov_matrix.shape[1]
else:
    cov_matrix = []
    num_covs = 0
print ("The num_covs is %r" % num_covs)
###########################
while not finished:
    try:
        test_snps = []
        # parse test SNP and associated info from input file row
        for i in range(len(infiles)):
            test_snps.append(test_parse_test_snp.parse_test_snp(snpinfo[i], options))
        rescale_totals.rescale_totals(test_snps)
        ref_as_counts = sum([np.sum(x.AS_target_ref) for x in test_snps])
        alt_as_counts = sum([np.sum(x.AS_target_alt) for x in test_snps])
        tot_as_counts = ref_as_counts + alt_as_counts
        all_counts = sum([test_snps[i].counts for i in range(len(test_snps))])
        if tot_as_counts < options.min_as_counts:
            # skip, not enough allele-specific counts
            for i in range(len(infiles)):
                line = infiles[i].readline().strip()
                if line:
                    snpinfo[i] = line.split()
                else:
                    finished = True
            continue
        
        row_count+=1
        old_genos = [test_snps[y].geno_hap1 + test_snps[y].geno_hap2
                     for y in range(len(test_snps))]
       
        if options.shuffle:
            # permute genotypes
            perm = range(len(test_snps))
            random.Random(options.seed).shuffle(perm)
            geno1temp = [test_snps[y].geno_hap1 for y in perm]
            geno2temp = [test_snps[y].geno_hap2 for y in perm]
            for i in range(len(test_snps)):
                test_snps[i].geno_hap1 = geno1temp[i]
                test_snps[i].geno_hap2 = geno2temp[i]
 
        starting_gene = [np.float64(x) for x in [0.1, 0.001]]
        maxlike = 10000000000
        for start in starting_gene:
            starts = [np.float64(options.initial), np.float64(options.initial), np.float64(0), np.float64(start)]
            res = minimize(ll.ll_one, starts, args=(disease_stat, test_snps, True, options.is_as_only, bnb_sigmas, as_sigmas,
                                                      options.read_error_rate, [], cov_matrix),
                           options={"maxiter" : 50000, "disp" : options.verbose},
                           method=OPTIMIZER)
            new_par = res.x
            new_loglike = ll.ll_one(new_par, disease_stat, test_snps, options.is_bnb_only, options.is_as_only, bnb_sigmas,
                                            as_sigmas, options.read_error_rate, [], cov_matrix)
            if new_loglike < maxlike:
               starting_par = new_par

        cov_coefs=[]
        for cov in range(num_covs):
            res = minimize(ll.ll_cov, [np.float64(0)],
                           args=(starting_par, disease_stat, test_snps, True, options.is_as_only, bnb_sigmas, as_sigmas,
                                 options.read_error_rate, cov_coefs, cov_matrix),
                           options={"maxiter" : 50000, "disp" : options.verbose},
                           method=OPTIMIZER)

            new_coef = res.x
            cov_coefs = np.concatenate([cov_coefs, new_coef])
        
        res = minimize(ll.ll_one, starting_par,
                       args=(disease_stat, test_snps, options.is_bnb_only, options.is_as_only, bnb_sigmas,
                             as_sigmas, options.read_error_rate, cov_coefs, cov_matrix),
                             options={"maxiter" : 50000, "disp" : options.verbose},
                             method=OPTIMIZER)
        best1par = res.x
        loglike1par = ll.ll_one(best1par, disease_stat, test_snps, options.is_bnb_only, options.is_as_only, bnb_sigmas,
                                        as_sigmas, options.read_error_rate, cov_coefs, cov_matrix) 
        start = [best1par[0], best1par[1], best1par[2], best1par[2], best1par[3]]
        
        res = minimize(ll.ll_two, start,
                       args=(disease_stat, test_snps, options.is_bnb_only, options.is_as_only, bnb_sigmas, as_sigmas, 
                             options.read_error_rate, cov_coefs, cov_matrix),
                             options={"maxiter" : 50000, "disp" : options.verbose},
                             method=OPTIMIZER)
        best2par = res.x
        loglike2par = ll.ll_two(best2par, disease_stat, test_snps, options.is_bnb_only, options.is_as_only, bnb_sigmas,
                                        as_sigmas, options.read_error_rate, cov_coefs, cov_matrix)
        
        write_results.write_results(outfile, snpinfo, loglike1par, loglike2par, best2par, tot_as_counts, ref_as_counts, alt_as_counts, all_counts)

    except Exception:
       # write_results.write_results(outfile, snpinfo, "-9", "-9", "-9", tot_as_counts, ref_as_counts, alt_as_counts, all_counts)
        write_empty_result(outfile, snpinfo)
        pass

    # read next set of lines from input file
    for i in range(len(infiles)):
        line = infiles[i].readline().strip()
        if line:
            snpinfo[i] = line.split()
        else:
            finished = True

