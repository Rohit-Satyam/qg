import glob
import os

index = '/home/parashar/scratch/hg19_resource/bowtie2_index/genome'

in_dirs = ['Na_PhenDC3/SRR12', 'Na_PhenDC3/SRR34']
# in_dirs = ['Na_K_1/SRR1693705', 'Na_K_1/SRR1693706']
# in_dirs = ['Na_K_2/SRR1693707', 'Na_K_2/SRR1693708']
# in_dirs = ['Na_PDS_1/SRR1693709', 'Na_PDS_1/SRR1693710']
# in_dirs = ['Na_PDS_2/SRR1693711', 'Na_PDS_2/SRR1693712']

out_dir = '../data/g4_seq/reanalysis/aligned_sequences'
src = 'bowtie2'
for d in in_dirs:
    od = out_dir + '/' + d
    fns = glob.glob('../data/g4_seq/reanalysis/raw_data/%s/control/*.fastq' % d)
    for fn in fns:
        f = d[-2:]+"_"+fn.split('/')[-1].split('.')[0]
        of = os.path.join(od, '%s.sam' % f)
        cmd = "bsub -q debugq -J %s -n 4 -R 'span[hosts=1]' -o ./bowtie_logs/%s.log '%s -p 1 -x %s -U %s -S %s'" % (f, f, src, index, fn, of)
        os.system(cmd)
        print (cmd)
