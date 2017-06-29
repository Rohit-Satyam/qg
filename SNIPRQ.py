import numpy as np
import networkx as nx
from collections import deque

# This code has been tested for Python 3.6 only.

class StemGroups(object):
    """
        Searches for putative stems of G4. Identifies all stems of
        pattern GGG, GNGG or GGNG starting from each base of the
        sequence, where length of N can be upto max_bulge_len. This
        function only considers stems with one bulge. Also only three
        Gs will be considered as quartet forming Gs.

        Input parameters:
            sequence: A nucleotide sequence (type: str)
            max_bulge_len: Maximum length of a bulge allowed in stem
                           (type: int)
            max_loop_len: Maximum length of a loop. i.e distance
                          between two stems (type: int)
            base: Can be either 'G' or 'C' depending which strand to
                  search (type: str)
    """

    def __init__(self, sequence, max_loop_len, max_bulge_len, base):

        self.stems = []
        self.stemLen = 3  # QuadGraph supports stem length of 3 only.
        self.sequence = sequence.upper()
        self.base = base
        self.maxBulgeLen = max_bulge_len
        self.maxLoopLen = max_loop_len
        self.kmerSize = self.stemLen + self.maxBulgeLen
        self.stems = []
        self.pos = None
        self.kmer = None
        self.k_stems = None

    def generate_stems(self):
        """
            A generators that yields groups of stems such that,
            within a group the stems are guaranteed to form a
            G4 motif and also stems from different groups will not
            form a G4 motif given the maximum allowed loop length

            Yield value:
                numpy array of shape s x 3, where s is the numnber of
                stems found. For each stem sequence position of each G
                within the stem is indicated.
        """

        for i in range(len(self.sequence) - self.kmerSize + 1):
            self.pos = i
            self.kmer = self.sequence[i: i + self.kmerSize]
            self.k_stems = self._get_stems_in_kmer(self.kmer, self.pos)
            if self._update_and_release_stems() is True:
                yield np.array(self.stems)
                self.stems = list(self.k_stems)

        # iterating over sub-kmers of last kmer
        for i in range(1, len(self.kmer) - self.stemLen + 1):
            self.k_stems = self._get_stems_in_kmer(self.kmer[i:],
                                             self.pos + i)
            if self._update_and_release_stems() is True:
                yield np.array(self.stems)
                self.stems = list(self.k_stems)

        if len(self.stems) >= 4:
            yield np.array(self.stems)

        self.stems = []

    def _get_stems_in_kmer(self, k, p):
        """
            Identifies the perfect stem if present or a single bulged
            stem without 'G' in the bulge.
            Returns:
                One stem in the kmer (type: list)
        """
        if k[:3] == self.base * 3:
            return [[p, p + 1, p + 2]]
        else:
            k = np.array(list(k))
            ks = len(k)
            if k[0] == self.base:
                if k[1] == self.base:
                    return [[p, p + 1, p + i] for i in
                            np.where(k == self.base)[0][2:3]]
                else:
                    for i in range(2, ks - 1):
                        if k[i] == self.base and k[i + 1] == self.base:
                            return [[p, p + i, p + i + 1]]
                            break
        return []

    def _update_and_release_stems(self):
        """
            Checks whether the set of stems for current kmer be
            merged to the current group of stems or the stems be
            yielded, the third scenario is where the stem group
            is reset as stem form the current kmers and previous
            stems are completely purged. This purging happens when
            either the stem group has less than four non overlappping
            stems and the new stems of current kmer are farther than
            the maximum allowed loop length.

        """
        len_stems = len(self.stems)
        if len(self.k_stems) > 0:
            if len_stems > 0:
                d = self.k_stems[0][0] - self.stems[-1][2]
                if d > self.maxLoopLen + 1:
                    if len_stems >= 4:
                        return True
                    else:
                        self.stems = list(self.k_stems)
                else:
                    self.stems.extend(self.k_stems)
            else:
                self.stems.extend(self.k_stems)
        return False


class QuadPaths(object):

    """
        Contructs a graph of stems such that all stems that are
        not more than user provided distance away from each other
        are connected by the edge. Thereafter, generates all
        possible 4 node paths.

        Input parameters:
            stems: A 2D array of stems and base positions (like the
                   one yielded by GetStems) (type: np.array)
            max_loop_len: Maximum allowed loop length (type: int)
    """

    def __init__(self, stems, max_loop_len):

        self.stems = stems
        self.maxLoopLen = max_loop_len
        self.graphStart = self.stems[0][0]
        self.nodeCodes = {}
        self.nodeInQuad = {}
        self.sortedNodeCodes = []
        self.G = self._populate_graph()
        self.QPS = np.array(self._make_quad_paths())
        self._remove_non_quad_nodes()
        if self.QPS.shape[0] > 0:
            self.info = [
                len(self.G.nodes()), len(self.G.edges()),
                self.QPS.shape[0], self.QPS[0][0][0], self.QPS[-1][-1][-1]
            ]
        else:
            self.info = []

    def _populate_graph(self):
        """
            Construct a directed acyclic graph of stems. Each edge
            represents the loop of a putative G4. Is assued that the
            stems are 2D array of shape nx3, where n is the number
            of stems and 3 represents the position of each G in the
            corresponding stem. The array should be pre-sorted by first
            stem position to ensure that the resultant graph is acyclic.
        """
        G = nx.DiGraph()
        all_starts = self.stems[:, 0]
        for stem in self.stems:
            targets = np.where((all_starts > stem[2] + 1) &
                     (all_starts <= stem[2] + self.maxLoopLen + 1))[0]
            n1 = self._stem_encoder(stem)
            self.nodeInQuad[n1] = False
            self.nodeCodes[n1] = stem
            self.sortedNodeCodes.append(n1)
            for t in targets:
                n2 = self._stem_encoder(self.stems[t])
                G.add_edge(n1, n2)
        return G

    def _stem_encoder(self, stem_array):
        """
            Represents the bases of a stem in a string notation.
            For example: '41*---*' represents a stem with three
            bulges and stems at position 41, 42 and 46.
        """
        diff = np.diff(stem_array)
        if diff[0] != 1:
            return '%d%s**' % (stem_array[0] - self.graphStart,
                               '-' * (diff[0] - 1))
        else:
            return '%d*%s*' % (stem_array[0] - self.graphStart,
                               '-' * (diff[1] - 1))

    def _make_quad_paths(self):
        """
            Recursive iteration over all 4-node paths in the
            graph to save all the quad-paths in the graph. This
            function also helps identify stems which do not form
            a quad-path.
        """
        quadpaths = deque()
        max_l = self.maxLoopLen + 1
        for n1 in self.sortedNodeCodes:
            if n1 not in self.G:
                continue
            n1_a = self.nodeCodes[n1]
            s2 = self.G.successors(n1)
            for n2 in self.sortedNodeCodes:
                if n2 not in s2:
                    continue
                n2_a = self.nodeCodes[n2]
                s3 = self.G.successors(n2)
                for n3 in self.sortedNodeCodes:
                    if n3 not in s3:
                        continue
                    n3_a = self.nodeCodes[n3]
                    s4 = self.G.successors(n3)
                    for n4 in self.sortedNodeCodes:
                        if n4 not in s4:
                            continue
                        n4_a = self.nodeCodes[n4]
                        for i in [n1, n2, n3, n4]:
                            self.nodeInQuad[i] = True
                        quadpaths.append((n1_a, n2_a, n3_a, n4_a))
        return quadpaths

    def _remove_non_quad_nodes(self):
        """
            Remove nodes that do not occur in any quad-path from
            the graph.
        """
        del_nodes = []
        for node in self.G.nodes():
            if self.nodeInQuad[node] is False:
                del_nodes.append(node)
        for node in del_nodes:
            self.G.remove_node(node)
        return True
