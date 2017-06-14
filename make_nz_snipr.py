
import numpy as np
import sys

chrom = sys.argv[1]
strand = sys.argv[2]

a = np.load('../data/snipr/%s_%s_scores.npy' % (chrom, strand))
np.save('../data/snipr/%s_%s_scores_NZ.npy' % (chrom, strand), a[a>0])