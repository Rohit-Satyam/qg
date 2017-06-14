
import os
import numpy as np
import glob

if __name__ == '__main__':
    chroms = ['chr'+str(x) for x in range(1,23)] + ['chrX', 'chrY']
    base_dir = "/home/parashar/scratch/quadcomb/data/REMAP_TFBS/tf_wise_files"
    tfs = np.array([x.split('/')[-1] for x in sorted(glob.glob("%s/*" % base_dir))])

    script = 'TFBS_snipr.py'
    for chrom in chroms:
        for tf in tfs:
            sig = "%s_%s" % (chrom, tf)
            cmd = 'bsub -q debugq -J %s -o ./tfbs_snipr_logs/%s.log python %s %s %s' % (
                    sig, sig, script, chrom, tf)
            print (cmd)
            os.system(cmd)