import glob
import os

if __name__ == "__main__":
    chroms = ['chr%d' % i for i in range(1, 23)] + ['chrX', 'chrY']
    strands = ["positive", 'negative']
    home_dir = "/home/parashar/scratch/quadcomb/data/quad_graphs_mll50"
    out_dir = '/home/parashar/scratch/quadcomb/data/quad_graph_loop_profile'
    script = 'QG_all_loop_len_calc.py'
    log_dir = 'qg_loops_log'
    n = 0
    for chrom in chroms:
        for strand in strands:
            base_dir = "%s/%s/%s" % (home_dir, chrom, strand)
            dirs = sorted(glob.glob("%s/*" % base_dir))
            for d in dirs:
                od = "%s/%s/%s" % (out_dir, chrom, strand)
                if not os.path.isdir(od):
                    os.makedirs(od)
                out_file = "%s/%s.json" % (od, d.split('/')[-1])
                sig = chrom[-2:] + strand[0] + d.split('/')[-1]
                bsub = "bsub -q debugq -o %s/%s -J %s" % (log_dir, sig, sig)
                cmd = "%s python %s %s %s" % (bsub, script, d, out_file)
                os.system(cmd)
                n += 1
print ("%d jobs submitted" % n)
