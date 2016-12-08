import glob
import os

if __name__ == "__main__":
    base_dir = '/home/parashar/scratch/quadcomb/data'
    for sample in ['Na_K_1', 'Na_K_2', 'Na_PDS_1', 'Na_PDS_2']:
        in_dir = "%s/QuadGraph_switchpoint_intersect/%s" % (base_dir, sample)
        out_dir = "%s/loop_lens/%s" % (base_dir, sample)
        bedfiles = glob.glob("%s/*" % in_dir)
        for bed_file in bedfiles:
            fn = bed_file.split('/')[-1]
            out_file = "%s/%s.json" % (out_dir, fn)
            bsub = "bsub -q debugq -J %s -o loop_lens_logs/%s.log" % (fn, fn)
            cmd = "%s python Loop_len_calc_g4seq_quadgraph.py %s %s" % (
                bsub, bed_file, out_file)
            os.system(cmd)
