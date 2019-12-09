def write_header(outfile):
    outfile.write("\t".join(["TEST.SNP.CHROM", "TEST.SNP.POS",
                             "LOGLIKE.NULL", "LOGLIKE.ALT",
                             "CHISQ", "P.VALUE", "ALPHA", "BETA",
                             "PHI", "TOTAL.AS.READ.COUNT",
                             "REF.AS.READ.COUNT", "ALT.AS.READ.COUNT",
                             "TOTAL.READ.COUNT"]) + "\n")

