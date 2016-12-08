import os

chrom_dir = '/home/parashar/scratch/hg19_resource/chromosomes'
#chroms = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY']
chroms = ['chr5', 'chr10']
out_dir = '/home/parashar/scratch/quadcomb/data/quad_graphs_mll50/'

for chrom in chroms:
    for base in ['G', 'C']:
        fasta = "%s/%s.fa" % (chrom_dir, chrom)
        bsub = "bsub -q shortq -J %s -n 16 -R 'span[hosts=1]' -o ./quad_graph_logs/%s.log" % (chrom, chrom)
        cmd = "%s python QuadGraph_parallel.py %s %s %s 3 50 15" % (bsub, fasta, out_dir, base)
        os.system(cmd)
