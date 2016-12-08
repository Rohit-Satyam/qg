import networkx as nx
import numpy as np
import sys
import json


class LoopFinder():
    def __init__(self, graphfile, pos, strand='+'):
        print (sorted((nx.read_graphml(graphfile).nodes())))
        self.stems = sorted(list(set([
            int(x.split('*')[0].split('-')[0]) for x in
            nx.read_graphml(graphfile).nodes() if x.count('-') == 0
        ])))
        self.pos = pos
        print (self.uniqueStemStarts)
        # self.ps = sorted([int(x.split('*')[0]) for x in
        #                   nx.read_graphml(graphfile).nodes() if
        #                   x.count('-') == 0])
        # self.ps = np.array([[x, x + 1, x + 2] for x in self.ps])

        # if self.ps.shape[0] > 0:
        if strand == '-':
            self.ps = self.ps[::-1, ::-1] * -1
            self.pos = -1 * self.pos

        self.pos_stems = [self._find_stem(x) for x in range(3)]
        self.next_pos_stems = [
            self._get_next_nov_stem(x) for x in self.pos_stems]
        self.stems_pos_loop_len = [
            self._get_loop_len(s1, s2) for s1, s2 in zip(
                self.pos_stems, self.next_pos_stems
            )]
        self.start_stem = self._get_next_stem()
        self.pre_start_pos_loop_len = self._get_loop_len(
            self.start_stem, self._get_next_nov_stem(self.start_stem))
        self.inter_stem_pos_loop_len = self._get_loop_len(
            self._get_prev_stem(), self.start_stem)

    def _find_stem(self, stempos):
        ret_val = self.ps[self.ps[:, stempos] == self.pos]
        if ret_val.shape[0] > 0:
            return list(ret_val[0])
        return []

    def _get_next_nov_stem(self, stem):
        if len(stem) != 3:
            return []
        ret_val = self.ps[self.ps[:, 0] > stem[2] + 1]
        if ret_val.shape[0] > 0:
            return list(ret_val[0])
        return []

    def _get_loop_len(self, s1, s2):
        if len(s1) != 3 or len(s2) != 3:
            return -1
        return s2[0] - s1[2] - 1

    def _get_next_stem(self):
        ret_val = self.ps[self.ps[:, 0] > self.pos]
        if ret_val.shape[0] > 0:
            return list(ret_val[0])
        return []

    def _get_prev_stem(self):
        ret_val = self.ps[self.ps[:, 2] < self.pos]
        if ret_val.shape[0] > 0:
            return list(ret_val[-1])
        return []

    def get_all_loop_lens(self):
        if self.ps.shape[0] > 0:
            return list(map(str, [
                self.stems_pos_loop_len[0],
                self.stems_pos_loop_len[1],
                self.stems_pos_loop_len[2],
                self.pre_start_pos_loop_len,
                self.inter_stem_pos_loop_len
            ]))
        else:
            return list(map(str, [-1, -1, -1, -1, -1]))


def get_graph_and_gss(bedfile):
    with open(bedfile) as h:
        bed = h.readlines()
    for bed_line in bed:
        c = str(bed_line).rstrip('\n').split('\t')
        start = int(c[3])
        base_dir = '/home/parashar/scratch/quadcomb/data/quad_graphs_mll50'
        strand = 'positive' if c[5] == '+' else 'negative'
        graph_dir = "%s/%s/%s/%s" % (
            base_dir, c[0], strand,
            int(((start // 1e6) * 1e6) + 1e6)
        )
        graph_file = "%s/%d.graphml" % (graph_dir, start)
        # switching strand because switchpoint and Quadgraph
        # were intersected on reverse strands
        strand = '+' if str(c[5]) == '-' else '-'
        yield graph_file, int(c[4]), strand


if __name__ == "__main__":
    bed_file = sys.argv[1]
    out_file = sys.argv[2]
    loops_length = []
    for n, i in enumerate(get_graph_and_gss(bed_file)):
        print("\r%d" % n, end='', flush=True)
        loops_length.append(LoopFinder(*i).get_all_loop_lens())
    with open(out_file, 'w') as OUT:
        json.dump(loops_length, OUT, indent=2)
