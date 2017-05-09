import os

if __name__ == "__main__":
    chroms = ['chr' + str(x) for x in range(1,23)] + ['chrX', 'chrY']
    log_dir = 'null_dist_PG4'
    script = 'PG4_null_dist_calc.py'
    n = 0
    for bulge in ['0', '5']:
        for chrom in chroms:
            for loop in ['5', '10', '15', '25']:
                for feat_len in ['50', '100', '500', '1000']:
                    sig = "%s_%s_%s_%s" % (chrom, loop, bulge, feat_len)
                    bsub = "bsub -q debugq -o %s/%s -J %s" % (log_dir, sig, sig)
                    cmd = "%s python %s %s %s %s %s" % (bsub, script, chrom, loop, bulge, feat_len)
                    print (cmd)
                    os.system(cmd)
                    n +=1
print ("%d jobs submitted" % n)