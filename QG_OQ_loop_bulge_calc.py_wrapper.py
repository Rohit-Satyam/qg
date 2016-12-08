import glob
import os

if __name__ == "__main__":
    n = 0
    base_dir = '/home/parashar/scratch/quadcomb/data'
    for sample in ['Na_K_1', 'Na_K_2', 'Na_PDS_1', 'Na_PDS_2']:
        in_dir = "%s/oq_switchpoints/%s" % (base_dir, sample)
        out_dir = "%s/QG_OQ_bulge_loop_min_lengths/%s" % (base_dir, sample)
        bedfiles = glob.glob("%s/*" % in_dir)
        for bed_file in bedfiles:
            fn = bed_file.split('/')[-1]
            out_file = "%s/%s.json" % (out_dir, fn)
            bsub = "bsub -q debugq -J %s -o bulge_loop_lens_log/%s.log" % (
                fn, fn)
            cmd = "%s python QG_OQ_loop_bulge_calc.py %s %s" % (
                bsub, bed_file, out_file)
            os.system(cmd)
            n += 1
print ("%d jobs submitted" % n)
