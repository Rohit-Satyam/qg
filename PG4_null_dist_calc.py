from tqdm import tqdm
import pybedtools as pbt
import numpy as np
import sys

if __name__ == '__main__':
    chrom = sys.argv[1]
    loop = int(sys.argv[2])
    bulge = int(sys.argv[3])
    feature_length = int(sys.argv[4])
    chrom_lens = {x.rstrip('\n').split('\t')[0]: int(x.rstrip('\n').split('\t')[1]) for x in
                  open('/home/parashar/scratch/hg19_resource/hg19.genome').readlines()}
    num_iters = 5000
    features_per_mb = 10
    
    save_name = '/home/parashar/scratch/quadcomb/data/null_dist_pg4s/%s_%d_%d_%d_%d' % (
                    chrom, loop, bulge, features_per_mb, feature_length)
    base_dir = '/home/parashar/scratch/quadruplexes/hg19'
    fn = '%s/g3_%d_%d_%s_nov.bed' % (base_dir, loop, bulge, chrom)
    g4_bed = pbt.BedTool(fn)
    num_regions = int(chrom_lens[chrom]/1e6)*features_per_mb
    
    all_counts = []
    for i in tqdm(range(num_iters)):
        starts = np.random.randint(0, chrom_lens[chrom]-feature_length, size=num_regions)
        ends = starts + feature_length
        feature_bed = pbt.BedTool('\n'.join(["\t".join([chrom, str(s), str(e)]) for s,e
                                         in zip(starts, ends)]), from_string=True)
        all_counts.append(feature_bed.intersect(g4_bed, u=True).count())
    with open(save_name, 'w') as OUT:
        OUT.write('\n'.join([str(x) for x in all_counts]))
