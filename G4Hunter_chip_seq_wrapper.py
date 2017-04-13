
import os

    
script = '/home/parashar/scratch/quadcomb/scripts/G4Hunter.r'
in_fn = '../data/chip_seq_g4/common_peaks.seq'
out_fn = '../data/chip_seq_g4/g4hunter_scores.txt'
bsub = "bsub -n 2 -q debugq -J peaks -o g4hunter_chipseq_logs/peaks.log"
cmd = '%s Rscript %s %s %s' % (bsub, script, in_fn, out_fn)
os.system(cmd)

for i in range(100):
    in_fn = '../data/chip_seq_g4/random_regions/random_peaks_%d.seq' % i
    out_fn = '../data/chip_seq_g4/g4hunter_random/random_%d.txt' % i
    bsub = "bsub -n 2 -q debugq -J %d -o g4hunter_chipseq_logs/%d.log" % (i, i)
    cmd = '%s Rscript %s %s %s' % (bsub, script, in_fn, out_fn)
    os.system(cmd)