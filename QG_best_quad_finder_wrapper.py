import glob
import os

if __name__ == "__main__":
    chroms = ['chr%d' % i for i in range(1, 23)] + ['chrX', 'chrY']
    strands = ["positive", 'negative']
    home_dir = "/home/parashar/scratch/quadcomb/data/quad_graphs_mll50"
    out_dir = '/home/parashar/scratch/quadcomb/data/quad_scores'
    script = 'QG_best_quad_finder.py'
    log_dir = 'quad_scores'
    n = 0
    for chrom in chroms:
        for strand in strands:
            base_dir = "%s/%s/%s" % (home_dir, chrom, strand)
            mega_dirs = sorted(glob.glob("%s/*" % base_dir))
            for d in mega_dirs:
                od = "%s/%s/%s" % (out_dir, chrom, strand)
                out_file = "%s/%s.json" % (od, d.split('/')[-1])
                if not os.path.isdir(od):
                    os.makedirs(od)
                sig = chrom[-2:] + strand[0] + d.split('/')[-1]
                bsub = "bsub -q debugq -o %s/%s -J %s" % (log_dir,
                                                          sig, sig)
                cmd = "%s python %s %s %s" % (bsub, script, d, out_file)
                os.system(cmd)
                print (n)
                n += 1

print ("\n%d jobs submitted" % n)
