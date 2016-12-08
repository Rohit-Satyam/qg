import networkx as nx
import sys
import json
from tqdm import tqdm


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


def get_oq_props(graphfile, switch_point_pos, strand):
    """
        Get minimum bulge and loop length at 0th or +1 swtichpoint
        from QuadGraph
    """
    if strand not in ['+', '-']:
        raise ValueError('Strand notation error')
    G = nx.read_graphml(graphfile)
    nodes = sorted(G.nodes())
    nodes_sp = []
    for node in nodes:
        n_start = int(node.split('-')[0].split('*')[0])
        if strand == "+":
            if n_start == (switch_point_pos + 1):
                nodes_sp.append(node)
        else:
            n_start = n_start + node.count('-') + 2
            if n_start == (switch_point_pos - 1):
                nodes_sp.append(node)
        if n_start == switch_point_pos:
            nodes_sp.append(node)
    min_bulge = 4
    min_loop = 51
    for node in nodes_sp:
        bulge = node.count('-')
        if bulge < min_bulge:
            min_bulge = bulge
        if bulge == 0:
            if strand == "+":
                edges = [x for x in G.out_edges_iter([node], data=True)]
            else:
                edges = [x for x in G.in_edges_iter([node], data=True)]
            for e in edges:
                check_len = False
                if strand == "+":
                    if e[1].count('-') == 0:
                        check_len = True
                else:
                    if e[0].count('-') == 0:
                        check_len = True
                if check_len:
                    if e[2]['length'] < min_loop:
                        min_loop = e[2]['length']
    if min_loop == 51:
        min_loop = None
    else:
        min_loop -= 1
    if len(nodes_sp) > 0:
        return (min_bulge, min_loop)
    return (None, None)


if __name__ == "__main__":
    bed_file = sys.argv[1]
    out_file = sys.argv[2]

    bulge_lengths = []
    loop_lengths = []
    for i in tqdm(get_graph_and_gss(bed_file)):
        res = get_oq_props(*i)
        if res[0] is not None:
            bulge_lengths.append(res[0])
            loop_lengths.append(res[1])
    with open(out_file, 'w') as OUT:
        json.dump([bulge_lengths, loop_lengths], OUT, indent=2)
