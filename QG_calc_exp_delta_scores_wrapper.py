import glob
import os

if __name__ == "__main__":
    chroms = ['chr%d' % i for i in range(1, 23)] + ['chrX', 'chrY']
    strands = ["positive", 'negative']
    base_dir = "/home/parashar/scratch/quadcomb/data"
    home_dir = "%s/quad_scores" % base_dir
    out_dir = '%s/quad_graph_single_base_scores' % base_dir
    script = 'QG_calc_exp_delta_scores.py'
    log_dir = 'base_scores'
    n = 0
    for chrom in chroms:
        for strand in strands:
            base_dir = "%s/%s/%s" % (home_dir, chrom, strand)
            in_files = sorted(glob.glob("%s/*" % base_dir))
            for fn in in_files:
                od = "%s/%s/%s" % (out_dir, chrom, strand)
                out_file = "%s/%s" % (od, fn.split('/')[-1].split('.json')[0])
                if not os.path.isdir(od):
                    os.makedirs(od)
                sig = chrom[-2:] + strand[0] + fn.split('/')[-1].split(
                    '.json')[0]
                bsub = "bsub -q debugq -o %s/%s -J %s" % (log_dir,
                                                          sig, sig)
                cmd = "%s python %s %s %s" % (bsub, script, fn, out_file)
                os.system(cmd)
                print (n)
                n += 1

print ("\n%d jobs submitted" % n)
