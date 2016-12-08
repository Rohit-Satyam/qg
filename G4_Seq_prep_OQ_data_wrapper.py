import os
import glob

sample_dirs = ['Na_PhenDC3/SRR12', 'Na_PhenDC3/SRR34']
# sample_dirs = ['Na_K_1/SRR1693705', 'Na_K_1/SRR1693706']
# sample_dirs = ['Na_K_2/SRR1693707', 'Na_K_2/SRR1693708']
# sample_dirs = ['Na_PDS_1/SRR1693709', 'Na_PDS_1/SRR1693710']
# sample_dirs = ['Na_PDS_2/SRR1693711', 'Na_PDS_2/SRR1693712']

base_dir = '/home/parashar/scratch/quadcomb/data/g4_seq/reanalysis'
log_dir = '/home/parashar/scratch/quadcomb/scripts/oq_data_logs'

n = 0
for sd in sample_dirs:
    control_fastq_dir = "%s/raw_data/%s/control" % (base_dir, sd)
    treated_fastq_dir = "%s/raw_data/%s/treated" % (base_dir, sd)
    aligned_dir = "%s/aligned_sequences/%s" % (base_dir, sd)
    out_dir_b = "%s/oq_bed/%s" % (base_dir, sd)
    out_dir_q = "%s/oq_qual_array/%s" % (base_dir, sd)
    out_dir_m = "%s/oq_mismatch_array/%s" % (base_dir, sd)
    for fn in glob.glob("%s/*.fastq" % control_fastq_dir):
        f = fn.split('/')[-1].split('.')[0]
        c = fn
        t = "%s/%s.fastq" % (treated_fastq_dir, f)
        a = "%s/%s_%s.sam" % (aligned_dir, fn.split('/')[-3][-2:], f)
        ob = "%s/%s.bed" % (out_dir_b, f)
        oq = "%s/%s" % (out_dir_q, f)
        om = "%s/%s" % (out_dir_m, f)
        bsub_line = "bsub -J %s_%s -o %s/%s_%s" % (
            fn.split('/')[-3][-2:], f, log_dir, fn.split('/')[-3][-2:], f)
        cmd = "%s python G4_Seq_prep_OQ_data.py %s %s %s %s %s %s" % (
            bsub_line, c, t, a, ob, oq, om)
        os.system(cmd)
        n += 1
print ('%d jobs submitted' % n)
