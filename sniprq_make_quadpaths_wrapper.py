import os

if __name__ == "__main__":

    chroms = ['chr'+str(x) for x in range(1,23)] + ['chrX', 'chrY', 'chrM']
    script = 'sniprq_make_quadpaths.py'
    logs_dir = 'snipr_make_qp_logs'
    for chrom in chroms:
        for base, strand in zip(['G', 'C'], ['positive', 'negative']):
            out_fn = "../data/quadpaths/%s_%s" % (chrom, strand)
            bsub = "bsub -q debugq -n 4 -R 'span[hosts=1]' -o %s/%s_%s.log" % (logs_dir, chrom, base)
            cmd = "%s python %s %s %s %s" % (bsub, script, chrom, base, out_fn)
            print (cmd)
            os.system(cmd)
