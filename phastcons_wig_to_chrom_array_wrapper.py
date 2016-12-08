import os

chroms = [str(x) for x in range(1, 23)] + ['X', 'Y']

for c in chroms:
	cmd = 'bsub -q shortq -n 16 -R "span[hosts=1]" -J %s -o phastcons_logs/%s.log python phastcons_wig_to_chrom_array.py %s' % (c, c, c)
	print (cmd)
