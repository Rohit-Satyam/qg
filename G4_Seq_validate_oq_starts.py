import sys
import pybedtools as pbt
import numpy as np
from tqdm import tqdm
from collections import Counter
from scipy import stats

cycles = 151 - 1
headcrop = 15 - 1


def get_validated_start(bed_file, switchpoint_file,
                        qual_file, mismatch_file):
    beds = pbt.BedTool(bed_file)
    switchpoints = np.load(switchpoint_file, mmap_mode='r')
    quals = np.load(qual_file, mmap_mode='r')
    mismatches = np.load(mismatch_file, mmap_mode='r')
    out_bed = []
    pre_qual = []
    post_qual = []
    for b, s, q, m in tqdm(zip(beds, switchpoints, quals, mismatches)):
        if np.nonzero(m) != len(m):
            max_p = 0
            pos = None
            for k, v in Counter(s[0]).items():
                if v > max_p:
                    pos = int(k)
                    max_p = v
            max_p = max_p / 2000
            q = [i if i > 0 else 0 for i in q]
            filter_level = 1
            if pos == 0 or pos == len(m):
                continue
            if 10 < pos < len(m) - 10:
                filter_level += 1
                if np.median(q[pos:]) - np.median(q[:pos]) > 10:
                    filter_level += 1
                    if stats.mannwhitneyu(q[pos:], q[:pos])[1] < 10e-5:
                        filter_level += 1
                        if sum(m[:pos]) < 2:
                            filter_level += 1
                            if max_p > 0.5:
                                filter_level += 1
                                adjust_factor = headcrop + pos
                                if b.strand == "+":
                                    start = b.start + adjust_factor
                                else:
                                    start = b.start + cycles - adjust_factor - 1 # -1 added later
                                out_bed.append("\t".join(map(str, [
                                    b.chrom, start, start + 1,
                                    b.name, 0, b.strand
                                ])))
            pre_qual.append([np.mean(q[:pos]), float(filter_level)])
            post_qual.append(np.mean(q[pos:]))
    return out_bed, pre_qual, post_qual


if __name__ == "__main__":
    bed_file = sys.argv[1]
    switchpoint_file = sys.argv[2]
    qual_file = sys.argv[3]
    mismatch_file = sys.argv[4]
    out_file = sys.argv[5]

    res = get_validated_start(bed_file, switchpoint_file,
                              qual_file, mismatch_file)
    out_bed, pre_qual, post_qual = res
    fn = out_file.split('.bed')[0]
    pbt.BedTool("\n".join(out_bed), from_string=True).sort().saveas(out_file)
    np.save("%s_pre_qual" % fn, np.array(pre_qual, dtype=float))
    np.save("%s_post_qual" % fn, np.array(post_qual))
