import os

chroms = list(reversed(['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY']))

for strand in ['positive', 'negative']:
    for chrom in chroms:
        bsub_line = "bsub -q debugq -R 'span[hosts=1]' -n 8"
        cmd = "%s python Graph_merge_jsons.py %s %s" % (bsub_line, chrom, strand)
        os.system(cmd)
