import numpy as np
import sys
from tqdm import tqdm

def bed_to_intervals(fn):
    intervals = []
    with open(fn) as h:
        for l in h:
            c = l.rstrip('\n').split('\t')
            intervals.append((int(c[1]), int(c[2])))
    return intervals

if __name__ == '__main__':
    
    chrom = sys.argv[1]
    tf = sys.argv[2]

    rand_scores = []
    tf_dir = "/home/parashar/scratch/quadcomb/data/REMAP_TFBS/tf_wise_files"
    snipr_dir = '/home/parashar/scratch/quadcomb/data/snipr'

    snipr = np.load('%s/%s_positive_scores.npy' % (snipr_dir, chrom), mmap_mode='r')
    snipr_nz = np.load('%s/%s_positive_scores_NZ.npy' % (snipr_dir, chrom), mmap_mode='r')
    nz_len = len(snipr_nz)
    
    scores = []
    rand_scores = []
    intervals = bed_to_intervals("%s/%s/%s.bed" % (tf_dir, tf, chrom))
    for interval in tqdm(intervals):
        s = snipr[interval[0]:interval[1]]
        s = s[s > 0]
        s_len = len(s)
        scores.append(s.sum())
        temp = []
        for i in range(1000):
            r = np.random.randint(0, nz_len-s_len)
            temp.append(snipr_nz[r:r+s_len].sum())
        rand_scores.append(temp)

    np.save('../data/REMAP_TFBS/snipr_tf/%s_%s_scores' % (chrom, tf), np.array(scores))
    np.save('../data/REMAP_TFBS/snipr_tf/%s_%s_random_scores' % (chrom, tf), np.array(rand_scores))
