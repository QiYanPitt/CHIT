import numpy as np
def read_as_sigmas(options, infiles):
    """Read overdispersion parameters for allele-specific test
    (Beta-Binomial). Expect one for each individual."""

    if (options.as_disp):
        disp_file = open(options.as_disp)
        line = disp_file.readline()
        as_sigmas = []
        while line:
            val = np.float64(line.strip())
            if val < 0.0 or val > 1.0:
                raise ValueError("expected as_sigma values to be "
                                 " in range 0.0-1.0, but got %g" %
                                 val)
            as_sigmas.append(np.float64(line.strip()))
            line = disp_file.readline()

        disp_file.close()

        if len(as_sigmas) != len(infiles):
            raise ValueError("expected %d values in as_disp file "
                             "(one for each input file) but got "
                             "%d" % (len(infiles), len(as_sigmas)))

    else:
        as_sigmas = [0.001] * len(infiles)

    return as_sigmas

