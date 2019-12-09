import numpy as np
def read_bnb_sigmas(options, infiles):
    """Read overdispersion parameters for beta-negative binomial.
    Expect one for each individual."""
    if (options.bnb_disp):
        disp_file = open(options.bnb_disp)
        line = disp_file.readline()
        bnb_sigmas = []
        while line:
            bnb_sigmas.append(np.float64(line.strip()))
            line = disp_file.readline()
        disp_file.close()

        if len(bnb_sigmas) != len(infiles):
            raise ValueError("expected %d values in bnb_disp file "
                             "(one for each input file) but got %d"
                             % (len(infiles), len(bnb_sigmas)))
    else:
        bnb_sigmas = [0.001]*len(infiles)

    return bnb_sigmas

