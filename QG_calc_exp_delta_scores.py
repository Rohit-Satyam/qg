import numpy as np
import json
import sys
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


def get_max_pos(quads, alts):
    quads = list(np.array(quads).flatten())
    vals = []
    for i in quads:
        vals.append(get_bases_from_stem(i)[2])
    alts = list(np.hstack(np.array(alts).flatten()))
    for i in alts:
        vals.append(get_bases_from_stem(i)[2])
    return max(vals)


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


if __name__ == "__main__":
    in_file = sys.argv[1]
    out_file_prefix = sys.argv[2]

    start_pos = int(in_file.split('/')[-1].split('.json')[0])
    data = json.load(open(in_file))

    exp_scores = []
    delta_scores = []

    for graph_pos in tqdm(data):
        quads, scores, alts = data[graph_pos]
        max_pos = get_max_pos(quads, alts)
        bscores = {x: [[], [], []] for x in range(0, max_pos + 1)}
        for quad, score, alt in zip(quads, scores, alts):
            for i in range(4):
                stem = quad[i]
                stem_bases = get_bases_from_stem(stem)
                for b in stem_bases:
                    bscores[b][0].append(score)
                alt_stems = alt[i]
                if len(alt_stems) > 1:
                    for alt_stem in alt_stems:
                        if alt_stem != stem:
                            alt_quad = sort_nodes(
                                quad[:i] + [alt_stem] + quad[i + 1:])
                            alt_score = calc_quad_score(alt_quad)
                            alt_bases = get_bases_from_stem(alt_quad[i])
                            delta_bases = set(stem_bases).difference(alt_bases)
                            for db in delta_bases:
                                bscores[db][1].append(alt_score)

        pos = int(graph_pos)
        for i in range(len(bscores)):
            if len(bscores[i][0]) > 0:
                exp_scores.append((pos + i, max(bscores[i][0])))
            if len(bscores[i][1]) > 0:
                delta_scores.append((pos + i, max(bscores[i][1])))

    np.save("%s_exp_score" % out_file_prefix, np.array(exp_scores))
    np.save("%s_delta_score" % out_file_prefix, np.array(delta_scores))
