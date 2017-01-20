import os
import numpy as np
import sys


def make_array(chrom):
    # hg19 coordinates
    chromInfo = {
        '1': 249250621, '2': 243199373, '3': 198022430,
        '4': 191154276, '5': 180915260, '6': 171115067,
        '7': 159138663, 'X': 155270560, '8': 146364022,
        '9': 141213431, '10': 135534747, '11': 135006516,
        '12': 133851895, '13': 115169878, '14': 107349540,
        '15': 102531392, '16': 90354753, '17': 81195210,
        '18': 78077248, '20': 63025520, 'Y': 59373566,
        '19': 59128983, '22': 51304566, '21': 48129895
    }
    return np.zeros(chromInfo[chrom])


def process_line(c):
    pos = int(c[1]) - 1
    if len(c[3]) == 1 and len(c[4]) == 1:
        info_line = c[7].split(';')
        for i in info_line:
            d = i.split('=')
            if d[0] == 'CAF':
                caf = float(d[1].split(',')[1])
                return (pos, caf)
    return False


def get_skip_lines_num(vcf_file):
    n = 0
    with open(vcf_file) as h:
        for l in h:
            if l[0] == '#':
                n += 1
            else:
                return n


def line_gen(vcf, chrom, skip):
    n = 0
    with open(vcf) as h:
        for i in range(skip):
            next(h)
        for l in h:
            c = l.rstrip('\n').split('\t')
            print ('\r%d' % n, end='')
            n += 1
            if c[0] == chrom:
                yield c


if  __name__ == '__main__':
    base_dir = '/home/parashar/scratch/quadcomb/data/dbsnp'
    vcf_file = "%s/downloaded_vcf/All_20160601.vcf" % base_dir
    out_dir = '%s/chrom_arrays' % base_dir
    chrom = sys.argv[1]

    chrom_array = make_array(chrom)
    skip_lines = get_skip_lines_num(vcf_file)

    for line in line_gen(vcf_file, chrom, skip_lines):
        res = process_line(line)
        if res is not False:
            chrom_array[res[0]] = res[1]

    print ('Saving array')
    np.save("%s/chr%s" % (out_dir, chrom), chrom_array)
