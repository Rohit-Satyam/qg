import os

if __name__ == '__main__':
    script = 'snipr_calculator_parallel.py'
    cmds = []
    for chrom in ['chr%d' % x for x in range(1, 3)] + ['chrX']:
        for base in ['G', 'C']:
            cmds.append('bsub -q smp -n 16 -R "span[hosts=1]" -o snipr_logs/%s_%s -J %s_%s python %s %s %s 15' % (
                    chrom, base, chrom, base, script, chrom, base))
    for chrom in ['chr%d' % x for x in range(3, 23)] + ['chrY', 'chrM']:
        for base in ['G', 'C']:
            cmds.append('bsub -q shortq -n 16 -R "span[hosts=1]" -o snipr_logs/%s_%s -J %s_%s python %s %s %s 15' % (
                    chrom, base, chrom, base, script, chrom, base))
    for cmd in cmds:
        os.system(cmd)

