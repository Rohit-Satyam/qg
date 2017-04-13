import numpy as np

chroms = ['chr'+str(x) for x in range(1,23)] + ['chrX', 'chrY']

for chrom in chroms:
    for strand in ['positive', 'negative']:
        print (chrom, strand)
        fn_in = '../data/QG_expectation_scores/%s_%s.npy' % (
            chrom, strand)
        fn_out = '../data/QG_expectation_scores_corrected_values/%s_%s' % (
            chrom, strand)
        a = np.load(fn_in)
        b = 2**(-1/a)
        np.save(fn_out, b)
        print (a.shape, b.shape)

