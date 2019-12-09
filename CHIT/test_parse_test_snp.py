import sys
import numpy as np
import random

class TestSNP:
    def __init__(self, name, geno_hap1, geno_hap2, AS_target_ref, AS_target_alt,
                 hetps, totals, counts):
        self.name = name
        self.geno_hap1 = geno_hap1
        self.geno_hap2 = geno_hap2
        self.AS_target_ref = AS_target_ref
        self.AS_target_alt = AS_target_alt
        self.hetps = hetps
        self.totals = totals
        self.counts = counts


    def is_het(self):
        """returns True if the test SNP is heterozygous"""
        return self.geno_hap1 != self.geno_hap2

    def is_homo_ref(self):
        """Returns True if test SNP is homozygous for reference allele"""
        return self.geno_hap1 == 0 and self.geno_hap2 == 0

    def is_homo_alt(self):
        """Returns True if test SNP is homozygous for non-reference allele"""
        return self.geno_hap1 == 1 and self.geno_hap2 == 1


def parse_test_snp(snpinfo, options):
    snp_id = snpinfo[2]
    if snpinfo[16] == "NA":
        # SNP is missing data
        tot = 0
    else:
        # these totals are later rescaled by dividing
        # by the minimum total across individuals to
        # put them into a reasonable range for
        # estimating alpha and beta
        tot = float(snpinfo[16])

    if snpinfo[6] == "NA":
        geno_hap1 = 0
        geno_hap2 = 0
    else:
        geno_hap1 = int(snpinfo[6].strip().split("|")[0])
        geno_hap2 = int(snpinfo[6].strip().split("|")[1])

    if snpinfo[15] == "NA":
        count = 0
    else:
        count = int(snpinfo[15])

    if snpinfo[9].strip() == "NA" or geno_hap1 == geno_hap2:
        # SNP is homozygous, so there is no AS info
        return TestSNP(snp_id, geno_hap1, geno_hap2, [], [], [], tot, count)
    else:
        # positions of target SNPs
        snp_locs = np.array([int(y.strip()) for y in snpinfo[9].split(';')])

        # counts of reads that match reference overlapping linked 'target' SNPs
        snp_as_ref = np.array([int(y) for y in snpinfo[12].split(';')])

        # counts of reads that match alternate allele
        snp_as_alt = np.array([int(y) for y in snpinfo[13].split(';')])

        # heterozygote probabilities
        snp_hetps = np.array([np.float64(y.strip())
                          for y in snpinfo[10].split(';')])

        # linkage probabilities, not currently used
        snp_linkageps = np.array([np.float64(y.strip())
                                  for y in snpinfo[11].split(';')])

        # same SNP should not be provided multiple times, this
        # can create problems with combined test. Warn and filter
        # duplicate SNPs
        uniq_loc, uniq_idx = np.unique(snp_locs, return_index=True)

        if options.dup_snp_warn and uniq_loc.shape[0] != snp_locs.shape[0]:
            sys.stderr.write("WARNING: discarding SNPs that are repeated "
                                     "multiple times in same line\n")
            options.dup_snp_warn = False

        snp_as_ref = snp_as_ref[uniq_idx]
        snp_as_alt = snp_as_alt[uniq_idx]
        snp_hetps = snp_hetps[uniq_idx]
        snp_linkageps = snp_linkageps[uniq_idx]

        if options.shuffle:
            # permute allele-specific read counts by flipping them randomly at
            # each SNP
            for y in range(len(snp_as_ref)):
                if random.randint(0, 1) == 1:
                    temp = snp_as_ref[y]
                    snp_as_ref[y] = snp_as_alt[y]
                    snp_as_alt[y] = temp

        return TestSNP(snp_id, geno_hap1, geno_hap2, snp_as_ref,
                       snp_as_alt, snp_hetps, tot, count)

