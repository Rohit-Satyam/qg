import glob
import os
import time


if __name__ == "__main__":
    base_dir = '/home/parashar/scratch/quadcomb/data/g4_seq/reanalysis'
    # dirs = [['Na_PhenDC3/SRR12', 'Na_PhenDC3/SRR34']]
    dirs = [['Na_K_1/SRR1693705', 'Na_K_1/SRR1693706'],
            ['Na_K_2/SRR1693707', 'Na_K_2/SRR1693708'],
            ['Na_PDS_1/SRR1693709', 'Na_PDS_1/SRR1693710'],
            ['Na_PDS_2/SRR1693711', 'Na_PDS_2/SRR1693712']]
    log_dir = "./validation_logs"
    script = "G4_Seq_validate_oq_starts.py"
    for sub_dirs in dirs:
        for sd in sub_dirs:
            bed_dir = base_dir + '/oq_bed/' + sd
            switchpoint_dir = base_dir + '/switchpoint_array/' + sd
            qual_dir = base_dir + '/oq_qual_array/' + sd
            mismatch_dir = base_dir + '/oq_mismatch_array/' + sd
            out_dir = base_dir + '/oq_validated_bed/' + sd
            bed_files = glob.glob("%s/*.bed" % (bed_dir))
            for b in bed_files:
                fn = b.split('/')[-1].split('.')[0]
                s = "%s/%s.npy" % (switchpoint_dir, fn)
                q = "%s/%s.npy" % (qual_dir, fn)
                m = "%s/%s.npy" % (mismatch_dir, fn)
                o = "%s/%s.bed" % (out_dir, fn)
                bsub = "bsub -n 1 -q debugq -J %s_%s -o %s/%s_%s " % (
                    sd[-2:], fn, log_dir, sd[-2:], fn)
                cmd = "%s python %s %s %s %s %s %s" % (
                    bsub, script, b, s, q, m, o)
                print (cmd)
                os.system(cmd)
                time.sleep(0.9)
