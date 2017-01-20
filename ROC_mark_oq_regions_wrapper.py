import glob
import os


samples = {
    'Na_K_1': ['SRR1693705', 'SRR1693706'],
    'Na_K_2': ['SRR1693707', 'SRR1693708'],
    'Na_PDS_1': ['SRR1693709', 'SRR1693710'],
    'Na_PDS_2': ['SRR1693711', 'SRR1693712'],
}
jobs = 0
base_in_dir = '/home/parashar/scratch/quadcomb/data/g4_seq/reanalysis'
out_base_dir = '/home/parashar/scratch/quadcomb/data/ROC_data/regions'
for sample in samples:
    for subdir in samples[sample]:
        sam_dir = '%s/aligned_sequences/%s/%s' % (
            base_in_dir, sample, subdir)
        sam_files = glob.glob("%s/*.sam" % sam_dir)
        for sf in sam_files:
            fn = sf.split('/')[-1].split('_')[-1].split('.')[0]
            control_fastq = '%s/raw_data/%s/%s/control/%s.fastq' % (
                base_in_dir, sample, subdir, fn)
            treated_fastq = '%s/raw_data/%s/%s/treated/%s.fastq' % (
                base_in_dir, sample, subdir, fn)
            out_file_prefix = '%s/%s/%s/%s' % (
                out_base_dir, sample, subdir, fn)
            sig = '%s_%s_%s' % (sample.split('_')[1] + sample.split('_')[2],
                                sf.split('/')[-1].split('_')[0], fn)
            bsub = "bsub -q debugq -J %s -o roc_regions_logs/%s" % (
                sig, sig)
            cmd = "%s python ROC_mark_oq_regions.py %s %s %s %s" % (
                bsub, sf, control_fastq, treated_fastq, out_file_prefix)
            print (cmd)
            os.system(cmd)
            jobs += 1

print ("Total %d jobs submitted" % jobs)
