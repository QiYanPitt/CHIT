def rescale_totals(test_snps):
    min_tot = min([s.totals for s in test_snps])

    if min_tot > 0:
        for s in test_snps:
            s.totals = s.totals / min_tot

