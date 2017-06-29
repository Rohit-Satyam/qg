import numpy as np
import sys


def get_read_names_from_sam(samfile, mapq_cutoff=30):
    names = {}
    with open(samfile) as h:
        while True:
            if next(h)[0:3] == '@PG':
                break
        for l in h:
            c = l.split('\t')
            if c[2] == "*":
                continue
            if int(c[4]) > mapq_cutoff:
                try:
                    ref = c[0].split('.')[1]  # <-
                except IndexError:
                    ref = c[0]
                strand = '+' if c[1] == '0' else '-'
                names[ref] = [c[2], c[3], strand]
    return names


def lazy_fastq_parser(f):
    with open(f) as h:
        for l in h:
            try:
                name = l.split(' ')[0].split('.')[1]
            except IndexError:
                name = l.split(' ')[0][1:]
            seq = next(h).rstrip('\n')
            next(h)
            qual = next(h).rstrip('\n')
            # only from 15th to till 99th cycle
            yield (name, seq[14:99], qual[14:99])


if __name__ == "__main__":
    control_fastq = sys.argv[1]
    treated_fastq = sys.argv[2]
    sam_file = sys.argv[3]
    bed_out_path = sys.argv[4]
    qual_out_path = sys.argv[5]
    mismatch_out_path = sys.argv[6]

    aligned_reads = get_read_names_from_sam(sam_file)
    fastq_iters = (
        lazy_fastq_parser(control_fastq),
        lazy_fastq_parser(treated_fastq),
    )

    qual_diffs = []
    mismatches = []
    oq_bed = []
    for c, t in zip(fastq_iters[0], fastq_iters[1]):
        if c[0] != t[0]:
            raise IndexError('Files are out of sync!')
        if c[0] in aligned_reads:
            m = [0 if i == j else 1 for i, j in zip(c[1], t[1])]
            if sum(m) > 30:
                mismatches.append(m)
                qual_diff = np.array([ord(x) - 33 for x in c[2]]) - \
                    np.array([ord(x) - 33 for x in t[2]])
                qual_diffs.append(qual_diff)
                temp = list(aligned_reads[c[0]])
                oq_bed.append("\t".join([
                    temp[0], str(int(temp[1]) - 1), temp[1],
                    c[0], '0', temp[2]
                ]))
    with open(bed_out_path, 'w') as OUT:
        OUT.write("\n".join(oq_bed))
    np.save(qual_out_path, np.array(qual_diffs, dtype=int))
    np.save(mismatch_out_path, np.array(mismatches, dtype=int))
    print ('DONE')
