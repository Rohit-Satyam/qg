import glob
import os
import re
from collections import Counter
import time

if __name__ == '__main__':
    base_dir = '/home/parashar/scratch/quadcomb/data/ENCODE_TFBS_Uniform'

    samples = {}
    with open("%s/bed_files_info.tsv" % base_dir) as h:
        for l in h:
            c = l.rstrip('\n').split('\t')
            attrs = {x.split('=')[0]: x.split('=')[1] for x in c[1].split('; ')}
            fn = c[0].split('.gz')[0]
            sample = "%s_%s" % (attrs['cell'], attrs['antibody'].split('_')[0])
            if attrs['quality'] == 'good':
                samples[sample] = fn
    
    for sample in samples:
        fn = "%s/bed_files/%s" % (base_dir, samples[sample])
        out_fn = "%s/exp_scores/%s" % (base_dir, sample)
        if os.path.isfile(fn):
            bsub = 'bsub -q debugq -J %s -o tfbs_logs/%s.log' % (sample, sample)
            cmd = "python TFBS_exp_scores.py %s %s" % (fn, out_fn)
            cmd = "%s %s"  % (bsub, cmd)
            os.system(cmd)
            time.sleep(20)
