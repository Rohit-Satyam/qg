
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.sandbox.stats.multicomp import multipletests
import sys

def load_scores(chrom, tf):
    a = np.load('../data/REMAP_TFBS/snipr_tf/%s_%s_scores.npy' % (chrom, tf),
                mmap_mode='r')
    b = np.load('../data/REMAP_TFBS/snipr_tf/%s_%s_random_scores.npy' % (chrom, tf),
               mmap_mode='r')
    return a,b

def calc_pval(a, b):
    pvals = []
    for i in b:
        r_a = [a[np.random.randint(len(a))] for x in range(1000)]
        r_b = [i[np.random.randint(len(i))] for x in range(1000)]
        pvals.append(mannwhitneyu(r_a, r_b)[1])       
    corrected_pvals = multipletests(pvals, alpha=0.05, method='holm',
                                    is_sorted=False, returnsorted=False)[1]
    return np.array([pvals, corrected_pvals])

if __name__ == '__main__':
    
    tf = sys.argv[1]
    
    scores = []
    rand_scores = []
    chroms = ['chr'+str(x) for x in range(1,23)] + ['chrX', 'chrY']
    for chrom in chroms:
        a, b = load_scores(chrom, tf)
        scores.extend(a)
        rand_scores.extend(b)
    scores = np.array(scores)
    rand_scores = np.array(rand_scores).T
    np.save('../data/REMAP_TFBS/pvals/%s.npy' % tf,
            calc_pval(scores, rand_scores))