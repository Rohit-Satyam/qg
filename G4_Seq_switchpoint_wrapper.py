import os
import glob

base_dir = '/home/parashar/scratch/quadcomb/data/g4_seq/reanalysis/'
n = 0

in_dirs = ['Na_PhenDC3/SRR12', 'Na_PhenDC3/SRR34']
# in_dirs = ['Na_K_1/SRR1693705', 'Na_K_1/SRR1693706']
# in_dirs = ['Na_K_2/SRR1693707', 'Na_K_2/SRR1693708']
# in_dirs = ['Na_PDS_1/SRR1693709', 'Na_PDS_1/SRR1693710']
# in_dirs = ['Na_PDS_2/SRR1693711', 'Na_PDS_2/SRR1693712']


for i in in_dirs:
    in_dir = base_dir + 'oq_mismatch_array/' + i
    out_dir = base_dir + 'switchpoint_array/' + i
    fns = glob.glob("%s/*.npy" % in_dir)
    for fn in fns:
        print (fn)
        #s = "%s_%s" % (in_dir[-2:], fn.split('/')[-1].split('_')[0])
        s = fn.split('/')[-1].split('.')[0]
        bsub_cmd = "bsub -n 1 -q debugq -J %s -o ./switchpoint_logs/%s.log" % (
            s, s)
        cmd = "%s python G4_Seq_switchpoint.py %s %s/%s" % (
            bsub_cmd, fn, out_dir, s)
        os.system(cmd)
        print (cmd)
        n += 1

print ("%d jobs submitted" % n)
