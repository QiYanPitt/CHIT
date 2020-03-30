import math
from scipy.optimize import *
from scipy.special import gammaln
from scipy.special import betaln
import scipy.stats
import numpy as np

# This is just ln(beta(a, b))
def betaln_asym(a,b):
    if b > a:
        a,b = b,a
    if a < 1e6:
        return betaln(a,b)
    l=gammaln(b)
    l -= b*math.log(a)
    l += b*(1-b)/(2*a)
    l += b*(1-b)*(1-2*b)/(12*a*a)
    l += -((b*(1-b))**2)/(12*a**3)
    return l

def BNB_loglike(k,mean,sigma,n):
    #Put variables in beta-NB form (n,a,b)
    #sys.stderr.write(str(sigma)+"\n")
    try:
        mean = max(mean,0.00001)
        logps = [math.log(n) - math.log(n + mean),
                 math.log(mean) - math.log(n + mean)]
    except:
        raise
        n_val=n
        pdb.set_trace()
    p=np.float64(n/(n+mean)) #This is the p in negative binomial

    if sigma < 0.00001: #> 18: #20:
        loglike=-betaln(n,k+1)-math.log(n+k)+n*logps[0]+k*logps[1]
        return loglike

    sigma=(1/sigma)**2 #+sigma*n
    sigma=sigma #+math.sqrt(sigma)/(p*(1-p))**2
    a = p*sigma+1
    b = (1-p)*sigma
    #Rising Pochhammer = gamma(k+n)/gamma(n)
    if k>0:
        loglike=-betaln_asym(n,k)-math.log(k)
    else:
        loglike=0
    #Add log(beta(a+n,b+k))
    loglike += betaln_asym(a+n,b+k)
    #Subtract log(beta(a,b))
    loglike -= betaln_asym(a,b)
    return loglike

def addlogs(loga, logb):
    """Helper function: perform numerically-stable addition in log space"""
    # When loga=3, logb=2 => a=exp(3), b=exp(2) => This is to calculate log(a+b)
    return max(loga, logb) + math.log(1 + math.exp(-abs(loga - logb)))

def AS_betabinom_loglike(logps, sigma, AS1, AS2, hetp, error):
    """Given parameter, returns log likelihood of allele-specific
    part of test. Note that some parts of equation have been
    canceled out"""
    a = math.exp(logps[0] + math.log(1/sigma**2 - 1))
    b = math.exp(logps[1] + math.log(1/sigma**2 - 1))
    part1 = 0
    part1 += betaln(AS1 + a, AS2 + b)
    part1 -= betaln(a, b)
    if hetp==1:
        return part1
    e1 = math.log(error) * AS1 + math.log(1 - error) * AS2
    e2 = math.log(error) * AS2 + math.log(1 - error) * AS1
    if hetp == 0:
        return addlogs(e1, e2)
    return addlogs(math.log(hetp)+part1, math.log(1-hetp) + addlogs(e1,e2))

def calc_cov_factor(cov_fits, covs, i):
    if len(cov_fits) > 0:
        return 1 + sum(cov_fits * covs[i,:len(cov_fits)])
    else:
        return 1

def loglikelihood(alpha, beta, par_disease_a, par_disease_b, r, disease_stat, test_snps, is_bnb_only,
                  is_as_only, bnb_sigmas, as_sigmas, error,
                  cov_coefs, cov_matrix):
    loglike = 0
    # if input values are outside of reasonable range return a
    # very high -loglike
    if alpha <= 0 or beta <= 0 or r <= 0 or r > 1:
        return 10000000
    
    for i in range(len(test_snps)):
        alpha2 = alpha*math.exp(par_disease_a*np.float64(disease_stat[i]))
        beta2 = beta*math.exp(par_disease_b*np.float64(disease_stat[i]))
        if(test_snps[i].is_homo_ref()):
            m = 2*alpha2*test_snps[i].totals * calc_cov_factor(cov_coefs, cov_matrix, i)
        elif(test_snps[i].is_homo_alt()):
            m = 2*beta2*test_snps[i].totals * calc_cov_factor(cov_coefs, cov_matrix, i)
        else:
            m = (alpha2+beta2)*test_snps[i].totals * calc_cov_factor(cov_coefs, cov_matrix, i)
        if m<0:
            m = 0.000001
        if not is_bnb_only:
            for j in range(len(test_snps[i].AS_target_ref)):
                if test_snps[i].hetps[j]>.9:
                    hetp = min(0.99, test_snps[i].hetps[j])
                    logps = [math.log(alpha2) - math.log(alpha2+beta2),
                             math.log(beta2) - math.log(alpha2+beta2)]
                    loglike += AS_betabinom_loglike(logps, as_sigmas[i],
                                                    test_snps[i].AS_target_ref[j],
                                                    test_snps[i].AS_target_alt[j],
                                                    hetp, error)
        if not is_as_only:
            l = BNB_loglike(test_snps[i].counts, m, r, bnb_sigmas[i]) # i is the the number infiles or subjects
            loglike += l
    return -loglike

def ll_one(x, disease_stat, test_snps, is_bnb_only, is_as_only, bnb_sigmas, as_sigmas, error, cov_coefs, cov_matrix):
    alpha = x[0]
    beta = x[1]
    par_disease_a = x[2]
    par_disease_b = x[2]
    r = x[3]
    return loglikelihood(alpha, beta, par_disease_a, par_disease_b, r, np.float64(disease_stat), test_snps, is_bnb_only,
                         is_as_only, bnb_sigmas, as_sigmas, error, cov_coefs, cov_matrix)

def ll_two(x, disease_stat, test_snps, is_bnb_only, is_as_only, bnb_sigmas, as_sigmas, error, cov_coefs, cov_matrix):
    alpha = x[0]
    beta = x[1]
    par_disease_a = x[2]
    par_disease_b = x[3]
    r = x[4]
    return loglikelihood(alpha, beta, par_disease_a, par_disease_b, r, np.float64(disease_stat), test_snps, is_bnb_only,
                         is_as_only, bnb_sigmas, as_sigmas, error, cov_coefs, cov_matrix)


def load_covariates(cov_file):
    infile=open(cov_file)
    cov_table=[]
    while True:
        line=infile.readline()
        if line:
            cov_table.append([np.float64(x) for x in line.strip().split()])
        else:
            break
    return np.array(cov_table, dtype=np.float64)

def ll_cov(x, params, disease_stat, test_snps, is_bnb_only, is_as_only, bnb_sigmas, as_sigmas, error, other_cov_coefs, cov_matrix):
    alpha = params[0]
    beta = params[1]
    par_disease_a = params[2]
    par_disease_b = params[2]
    r = params[3]
    cov_coefs=np.concatenate([other_cov_coefs,x])
    return loglikelihood(alpha, beta, par_disease_a, par_disease_b, r, np.float64(disease_stat), test_snps, is_bnb_only,
                         is_as_only, bnb_sigmas, as_sigmas, error, cov_coefs,cov_matrix)

