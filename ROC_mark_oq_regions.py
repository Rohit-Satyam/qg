import sys
import os


def get_aligned_reads(sam_file):
    print ("reading SAM file", flush=True)
    aligned_reads = {}
    with open(sam_file) as h:
        for l in h:
            if l[0] == '@':
                continue
            else:
                c = l.split('\t')
                if c[4] == '42':
                    aligned_reads[c[0]] = (c[2], int(c[3]) - 1, c[1])
    return aligned_reads


def get_regions(control_fastq, treated_fastq, aligned_reads):
    print ("Fetching regions", flush=True)
    chroms = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY']
    oq_regions = {}
    no_oq_regions = {}
    for chrom in chroms:
        oq_regions[chrom] = {}
        no_oq_regions[chrom] = {}
        for strand in ['+', '-']:
            oq_regions[chrom][strand] = []
            no_oq_regions[chrom][strand] = []

    hc = open(control_fastq)
    ht = open(treated_fastq)

    for lc, lt in zip(hc, ht):
        read_name = lc.split(' ')[0][1:]
        if read_name not in aligned_reads:
            for i in range(3):
                next(hc)
                next(ht)
            continue
        chrom = aligned_reads[read_name][0]
        if chrom == 'chrM':
            for i in range(3):
                next(hc)
                next(ht)
            continue
        start_pos = aligned_reads[read_name][1]
        strand = '+' if aligned_reads[
            read_name][2] == '0' else '-'
        seq_c = next(hc).rstrip('\n')
        seq_t = next(ht).rstrip('\n')
        positions = []
        mm = [1 if seq_c[x] != seq_t[x] else 0 for x in range(14, 99)]
        if strand == '+':
            positions = [start_pos + x for x in range(14, 99)]
        else:
            positions = [start_pos + 150 - x for x in range(14, 99)][::-1]
        if sum(mm) > 30:
            s = mm.index(1)
            fp = positions[s]
            lp = positions[-1] + 1
            oq_regions[chrom][strand].append([fp, lp])
        elif sum(mm) < 3:
            fp = positions[0]
            lp = positions[-1] + 1
            no_oq_regions[chrom][strand].append([fp, lp])
        for i in range(2):
            next(hc)
            next(ht)
    hc.close()
    ht.close()

    return oq_regions, no_oq_regions


def make_bed(regions, prefix):
    print ("making BED format file: %s" % prefix, flush=True)
    chroms = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY']
    bed = []
    for chrom in chroms:
        for strand in ['+', '-']:
            for i in regions[chrom][strand]:
                bed.append('\t'.join([
                    chrom, str(i[0]), str(i[1]),
                    '.', '1', strand
                ]))
    temp_file = "%s_temp.bed" % prefix
    sorted_file = "%s_sorted.bed" % prefix
    merged_file = "%s.bed" % prefix
    with open(temp_file, 'w') as OUT:
        OUT.write('\n'.join(bed))
    os.system("sort -k1,1 -k2,2n %s > %s" % (temp_file, sorted_file))
    os.system("bedtools merge -i %s -c 5 -o sum -s > %s" % (
        sorted_file, merged_file))
    os.system("rm %s %s" % (temp_file, sorted_file))
    return True


if __name__ == '__main__':
    sam_file = sys.argv[1]
    control_fastq = sys.argv[2]
    treated_fastq = sys.argv[3]
    out_file_prefix = sys.argv[4]

    oq_regions, no_oq_regions = get_regions(
        control_fastq, treated_fastq, get_aligned_reads(sam_file))
    make_bed(oq_regions, "%s_oq_regions" % out_file_prefix)
    make_bed(no_oq_regions, "%s_no_oq_regions" % out_file_prefix)
