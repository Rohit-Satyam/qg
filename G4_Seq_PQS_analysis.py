import pybedtools as pbt
import re
import sys
from tqdm import tqdm
import numpy as np


class QuadMotifFinder():

    def __init__(self, sequence, stem=3, loop_start=1,
                 loop_stop=7, greedy=True, bulge=0, strand='-'):
        self.sequence = sequence
        self.stemSize = stem
        self.minLoopSize = loop_start
        self.maxLoopSize = loop_stop
        self.greedySearch = greedy
        self.bulgeSize = bulge
        self.strand = strand
        self._sanitize_params()
        signature, motifOv, motifNov = self._make_motif()
        self.resOv = self._run_regex(motifOv, signature)
        self.resNov = self._run_regex(motifNov, signature)

    def _sanitize_params(self):
        try:
            self.stemSize = int(self.stemSize)
            self.minLoopSize = int(self.minLoopSize)
            self.maxLoopSize = int(self.maxLoopSize)
            self.bulgeSize = int(self.bulgeSize)
        except:
            raise TypeError(
                'Please provide only integer values for size paramaters')
        try:
            self.greedySearch = bool(self.greedySearch)
        except:
            raise TypeError(
                'Please provide only boolean values for greedySearch')
        if self.minLoopSize > self.maxLoopSize:
            raise ValueError(
                'Minimum loop size has to be smaller than maximum loop size')
        if self.bulgeSize > 0 and self.stemSize != 3:
            print ("Bulge search disabled as stem size not equal to 3")
            self.bulgeSize = 0

    def _make_base(self, strand):
        if strand == "+":
            base = 'G'
            invert_base = 'C'
        else:
            base = 'C'
            invert_base = 'G'
        return base, invert_base

    def _make_signature(self, base):
        return "".join(map(str, [base, self.stemSize, 'L',
                                 self.minLoopSize, '-', self.maxLoopSize]))

    def _make_stem_motif(self, base, invert_base):
        if self.bulgeSize > 0:
            stem_motif = '(?:%s|%s[AT%sN]{1,%d}%s|%s[AT%sN]{1,%d}%s)' % (
                base * 3, base, invert_base, self.bulgeSize, base * 2,
                base * 2, invert_base, self.bulgeSize, base)
        else:
            stem_motif = '%s' % (base * self.stemSize)
        return stem_motif

    def _make_loop_motif(self, base, invert_base):
        if self.greedySearch is True:
            loop_motif = "[ACTGN]{%d,%d}" % (
                self.minLoopSize, self.maxLoopSize)
        else:
            loop_motif = '[ATGCN](?:[ATN%s]|(?!%s{2})%s){%d,%d}' % (
                invert_base, base, base,
                self.minLoopSize - 1, self.maxLoopSize - 1)
        return loop_motif

    def _make_motif(self):
        base, invert_base = self._make_base(self.strand)
        quad_signature = self._make_signature(base)
        stem_motif = self._make_stem_motif(base, invert_base)
        loop_motif = self._make_loop_motif(base, invert_base)
        quad_motif_nov = (stem_motif + loop_motif) * 3 + stem_motif
        quad_motif_nov = "(" + quad_motif_nov + ")"
        quad_motif_ov = "(?=(%s))" % quad_motif_nov
        return quad_signature, quad_motif_ov, quad_motif_nov

    def _run_regex(self, pattern, signature):
        regex = re.compile(pattern, flags=re.IGNORECASE)
        reg_obj = regex.finditer(self.sequence)
        res = []
        for match in reg_obj:
            temp_res = [
                'QB', match.start(), match.start() + len(match.group(1)) + 1,
                len(match.group(1)), signature, match.group(1)
            ]
            res.append("\t".join(map(str, temp_res)))
        return res


if __name__ == "__main__":
    switch_point_bed = sys.argv[1]
    out_file = sys.argv[2]
    loop_start = int(sys.argv[3])
    loop_stop = int(sys.argv[4])

    offset = 50
    chrom_info = '/home/parashar/scratch/hg19_resource/hg19.genome'
    chrom_loc = '/home/parashar/scratch/hg19_resource/chromosomes'
    distances = []
    sequences = []
    for bedline in tqdm(pbt.BedTool(switch_point_bed).sort()):
        c = str(bedline).rstrip('\n').split('\t')
        if c[0] == 'chrM':
            continue
        bed = pbt.BedTool("\t".join(c), from_string=True).slop(
            b=offset, g=chrom_info)
        fasta_file = '%s/%s.fa' % (chrom_loc, c[0])
        bed = bed.sequence(fi=fasta_file, s=True)
        seq = [x.rstrip('\n').upper() for x in
               open(bed.seqfn).readlines()[1:]][0]
        q = QuadMotifFinder(seq, loop_start=loop_start,
                            loop_stop=loop_stop).resNov
        if len(q) > 0:
            q = [x.split('\t') for x in q][0]
            distances.append(int(q[1]))
            sequences.append(seq)
    np.save(out_file, np.array(distances))
    with open("%s.seq" % out_file, 'w') as OUT:
        OUT.write("\n".join(sequences))
