import numpy as np
import networkx as nx
import json
import sys
import glob
from tqdm import tqdm


def get_bases_from_stem(stem):
    us = int(stem.split('-')[0].split('*')[0])
    g1 = us
    ul = stem.count('-')
    if ul == 0:
        g2 = us + 1
        g3 = us + 2
    elif stem.split('-')[-1] == '**':
        g2 = us + ul + 1
        g3 = us + ul + 2
    else:
        g2 = us + 1
        g3 = us + ul + 2
    return [g1, g2, g3]


def sort_nodes(nodes):
    if nodes == []:
        return []
    a = []
    for x in nodes:
        a.append(get_bases_from_stem(x))
    a = np.array(a, dtype=int)
    sorted_index = np.lexsort((a[:, 2], a[:, 1], a[:, 0]))
    return [nodes[x] for x in sorted_index]


def calc_quad_score(quad):

    def calc_loop_size(u, v):
        us = int(u.split('-')[0].split('*')[0])
        vs = int(v.split('-')[0].split('*')[0])
        return vs - (us + 2) - 1

    def calc_loop_strength(loop_len, decay_rate=1.02):
        return np.exp(-decay_rate * loop_len)

    def calc_stem_strength(stem, decay_rate=1.9):
        bulge_size = stem.count('-')
        return np.exp(-decay_rate * bulge_size)

    def calc_gmean(vals):
        prod = 1
        for i in vals:
            prod = prod * i
        return prod**(1 / len(vals))

    def calc_hmean(l):
        return len(l) / sum([1 / x for x in l])

    loop_sizes = []
    stem_strengths = []
    for i in range(0, 3):
        loop_sizes.append(calc_loop_size(quad[i], quad[i + 1]))
    loop_strengths = [calc_loop_strength(ll) for ll in loop_sizes]
    stem_strengths = [calc_stem_strength(s) for s in quad]
    raw_scores = calc_gmean(
        [calc_gmean(loop_strengths), calc_hmean(stem_strengths)])
    transformed_scores = 1 / (-np.log2(raw_scores))
    return transformed_scores


def get_best_quads(graph):
    quads = []
    scores = []
    best_triplet = []
    for n1 in sort_nodes(graph.nodes()):
        if best_triplet != []:
            n1_stop = int(n1.split('*')[0].split('-')[0]) + 2 + n1.count('-')
            if int(best_triplet[0].split('*')[0].split('-')[0]) > n1_stop + 1:
                quad = [n1] + best_triplet
                quads.append(quad)
                scores.append(calc_quad_score(quad))
                continue
        best_score = 0
        best_quad = []
        for n2 in sort_nodes(graph.successors(n1)):
            for n3 in sort_nodes(graph.successors(n2)):
                for n4 in sort_nodes(graph.successors(n3)):
                    quad = [n1, n2, n3, n4]
                    score = calc_quad_score(quad)
                    if score > best_score:
                        best_quad = quad
                        best_score = score
        if best_quad != []:
            quads.append(best_quad)
            scores.append(best_score)
            best_triplet = best_quad[1:]
    return quads, scores


def get_quad_alternatives(graph, qs):
    alts = []
    for q in qs:
        alt = []
        alt.append(sort_nodes(graph.predecessors(q[1]) +
                              graph.successors(q[3])))
        alt.append(sort_nodes(list(set(
            graph.successors(q[0])).intersection(graph.predecessors(q[2])))))
        alt.append(sort_nodes(list(set(
            graph.successors(q[1])).intersection(graph.predecessors(q[3])))))
        alt.append(sort_nodes(graph.predecessors(q[0]) +
                              graph.successors(q[2])))
        alts.append(list(alt))
    return alts


if __name__ == "__main__":
    in_dir = sys.argv[1].rstrip('/')
    out_file = sys.argv[2]

    fns = sorted(glob.glob("%s/*.graphml" % (in_dir)))
    out_data = {}
    for fn in tqdm(fns):
        G = nx.read_graphml(fn)
        opt_quads, opt_scores = get_best_quads(G)
        alt_stems = get_quad_alternatives(G, opt_quads)

        start_pos = fn.split('/')[-1].split('.graphml')[0]
        out_data[start_pos] = [
            opt_quads,
            opt_scores,
            alt_stems,
        ]
    with open(out_file, 'w') as OUT:
        json.dump(out_data, OUT, indent=2)
