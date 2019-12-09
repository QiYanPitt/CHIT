import argparse

def parse_options():
    parser=argparse.ArgumentParser()
    parser.add_argument("-a", "--as_only",
                        action='store_true',
                        dest='is_as_only', default=False,
                        help="only perform the allele-specific part (Beta Binomial) "
                        "part of the test")
    #The store_true option automatically creates a default value of False. 
    parser.add_argument("-d", "--bnb_only", action='store_true',
                        dest='is_bnb_only', default=False,
                        help="only perform the association (Beta Negative Binomial) part "
                        "of the test")

    parser.add_argument("--cov_file", action='store',
                        dest='cov_file',
                        help="file containing covariates to include in the model"
                        ,default=None)

    parser.add_argument("-b", "--bnb_disp", action='store', dest='bnb_disp',
                        help="file containing depth (Beta Negative Binomial)"
                        "dispersion parameters", default=None)

    parser.add_argument("-o", "--as_disp", action='store',
                        dest='as_disp',
                        help="file containing allele-specific (Beta Binomial) dispersion "
                        "parameters", default=None)

    parser.add_argument("-s", "--shuffle", action='store_true',
                        dest='shuffle', default=False,
                        help="permute genotypes")

    parser.add_argument("-e", "--read_error_rate", action='store', dest='read_error_rate',
                        help="estimate of error rate, used to update "
                        "heterozygous genotype probabilities "
                        "(currently this option disabled / not used)",
                        type=float, default=0.005)

    parser.add_argument("-m", "--min_as_counts", action='store', dest='min_as_counts',
                        type=int, default=0,
                        help="only perform test when total number of allele-specific "
                        "read counts across individuals > MIN_COUNTS")

    parser.add_argument("-v", "--verbose", action='store_true', dest='verbose',
                        default=False, help="print extra information")

    parser.add_argument("--benchmark", dest="benchmark",
                        help="write information about time test is takes, number of optimization "
                        "functions, etc. to specified filename, or to stderr if '-' is specified")

    parser.add_argument("-q", "--equal", action='store_true',
                        dest='equal', default=False,
                        help="assume eaqual phenotype effect on both haplotypes")

    parser.add_argument("-i", "--initial", action='store', dest='initial',
                        type=float, default=0.5,
                        help="initial value for parameters alpha and beta")

    parser.add_argument("-r", "--seed", action='store', dest='seed',
                        type=int, default=1,
                        help="seed for random shuffle")

    parser.add_argument("infile_list", action='store', default=None)
    parser.add_argument("out_file", action='store', default=None)

    return parser.parse_args()

