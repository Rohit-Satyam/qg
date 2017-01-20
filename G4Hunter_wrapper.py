import glob
import os

base_dir = '/home/parashar/scratch/quadcomb/data/ROC_data'
for in_file in glob.glob("%s/seq_g4_hunter/*.seq" % base_dir):
    fn = in_file.split('/')[-1].split('.seq')[0]
    out_file = "%s/scores_g4_hunter/%s" % (base_dir, fn)
    bsub = 'bsub -n 8 -q debugq -J %s -o G4Hunter_logs/%s.log' % (fn, fn)
    cmd = "%s Rscript G4Hunter.r %s %s" % (bsub, in_file, out_file)
    os.system(cmd)
