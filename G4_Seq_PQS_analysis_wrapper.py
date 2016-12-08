import glob
import os
import time


if __name__ == "__main__":
    samples = {
        'Na_K_1': ['SRR1693705', 'SRR1693706'],
        'Na_K_2': ['SRR1693707', 'SRR1693708'],
        'Na_PDS_1': ['SRR1693709', 'SRR1693710'],
        'Na_PDS_2': ['SRR1693711', 'SRR1693712'],
    }
    base_dir = '/home/parashar/scratch/quadcomb/data/g4_seq/reanalysis'
    config = [
        (1, 3), (1, 7), (1, 15)
    ]
    script_name = 'G4_Seq_PQS_analysis.py'
    n = 0
    for sample in samples:
        for sub_dir in samples[sample]:
            switch_point_dir = "%s/oq_validated_bed/%s/%s" % (
                base_dir, sample, sub_dir)
            out_dir = "%s/oq_pq_distances/%s/%s" % (
                base_dir, sample, sub_dir)
            sp_bedfile = glob.glob("%s/*.bed" % switch_point_dir)
            for bf in sp_bedfile:
                fn = bf.split('/')[-1].split('.bed')[0]
                for c in config:
                    suffix = "_".join(list(map(str, c)))
                    out_file = "%s/%s_%s" % (out_dir, fn, suffix)
                    bsub = "bsub -q debugq -J %s_%s_%s -o %s/%s_%s_%s" % (
                        sub_dir[-2:], fn, suffix,
                        'oq_pq_distance_log', sub_dir[-2:], fn, suffix)
                    cmd = "%s python %s %s %s %d %d" % (
                        bsub, script_name, bf, out_file, c[0], c[1])
                    n += 1
                    print (cmd)
                    os.system(cmd)
                    time.sleep(0.005)
print ("%d jobs submitted" % n)
