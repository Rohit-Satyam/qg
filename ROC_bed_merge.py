import glob
import os
from tqdm import tqdm


def merge(in_fns, out_name):
    chroms = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY']
    out_bed = {}
    for chrom in chroms:
        out_bed[chrom] = []
    for fn in tqdm(in_fns, desc='Aggregating'):
        with open(fn) as h:
            for l in h:
                c = l.rstrip('\n').split('\t')
                out_bed[c[0]].append("\t".join([
                    c[0], c[1], c[2], '.', c[4], c[3]
                ]))
        h.close()

    for chrom in tqdm(chroms, desc='Merging'):
        temp_file = "%s_%s_temp.bed" % (out_name, chrom)
        with open(temp_file, 'w') as OUT:
            OUT.write("\n".join(out_bed[chrom]))

        sorted_file = "%s_%s_sorted.bed" % (out_name, chrom)
        os.system("sort -k2,2n %s > %s" % (temp_file, sorted_file))

        merged_file = "%s_%s.bed" % (out_name, chrom)
        os.system("bedtools merge -i %s -c 5 -o sum -s > %s" % (
            sorted_file, merged_file))
        os.system("rm %s %s" % (temp_file, sorted_file))
    return True


if __name__ == "__main__":
    samples = {
        'Na_K_1': ['SRR1693705', 'SRR1693706'],
        'Na_K_2': ['SRR1693707', 'SRR1693708'],
        'Na_PDS_1': ['SRR1693709', 'SRR1693710'],
        'Na_PDS_2': ['SRR1693711', 'SRR1693712'],
    }

    base_dir = '/home/parashar/scratch/quadcomb/data/ROC_data/regions'

    for sample in samples:
        print (sample)

        no_oq_fns = glob.glob("%s/%s/*/*_no_oq_regions.bed" % (
            base_dir, sample))
        no_oq_out_bed = "%s/chrom_wise/%s_no_oq_regions" % (base_dir, sample)
        merge(no_oq_fns, no_oq_out_bed)
        oq_fns = [x for x in glob.glob("%s/%s/*/*_oq_regions.bed" % (
            base_dir, sample)) if x not in no_oq_fns]
        oq_out_bed = "%s/chrom_wise/%s_oq_regions" % (base_dir, sample)
        merge(oq_fns, oq_out_bed)
