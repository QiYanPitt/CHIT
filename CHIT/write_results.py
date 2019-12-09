import scipy.stats

def write_results(outfile, snpinfo, loglike1par, loglike2par,
                  best2par, tot_as_counts, ref_as_counts, alt_as_counts,
                  all_counts):
    """Write result to output file. Tab-delimited columns are:
      1. chromosome,
      2. SNP position,
      3. Log likelihood 1 parameter model (Null)
      4. Log likelihood 2 parameter model (Alternative)
      3. Chi-squared statistic,
      4. P-value
      5. alpha parameter estimate (expression level
         of reference allele)
      6. beta parameter estimate (expression level of
         alternative allele)
      7. phi parameter estimate (beta-negative-binomial
         overdispersion
         parameter for this region)
      8. total number of allele-specific read counts for this
         region summed across individuals
      9. total number of reference haplotype allele-specific read counts
     10. total number of alt haplotype allele-specific read counts
     11. total number of mapped reads for this region,
         summed across individuals"""

    # compute likelihood ratio test statistic:
    chisq = 2 * (loglike1par - loglike2par)
    pval = (1-scipy.stats.chi2.cdf(chisq,1)),

    outfile.write("\t".join([snpinfo[0][0], snpinfo[0][1],
                             "%.2f" % -loglike1par,
                             "%.2f" % -loglike2par,
                             "%.3f" % chisq,
                             "%g" % pval,
                             "%g" % best2par[0],
                             "%g" % best2par[1],
                             "%g" % best2par[2],
                             "%d" % tot_as_counts,
                             "%d" % ref_as_counts,
                             "%d" % alt_as_counts,
                             "%d" % all_counts]) + '\n')
    outfile.flush()

