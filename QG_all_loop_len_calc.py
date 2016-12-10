import sys
from glob import glob
from tqdm import tqdm
from networkx import read_graphml
import json


if __name__ == "__main__":
    in_dir = sys.argv[1].rstrip('/')
    out_file = sys.argv[2]
    keys = [
        "0_0", "0_1", "0_2", "0_3",
        "1_0", "1_1", "1_2", "1_3",
        "2_0", "2_1", "2_2", "2_3",
        "3_0", "3_1", "3_2", "3_3",
    ]
    data = {}
    for k in keys:
        data[k] = {x: 0 for x in range(1, 51)}
    bulge_data = {0: 0, 1: 0, 2: 0, 3: 0}
    fns = glob("%s/*.graphml" % in_dir)
    for fn in tqdm(fns):
        G = read_graphml(fn)
        for i in G.edges_iter(data=True):
            k = "%d_%d" % (i[0].count('-'), i[1].count('-'))
            data[k][i[2]['length'] - 1] += 1
        for i in G.nodes():
            bulge_data[i.count('-')] += 1

    with open(out_file, 'w') as OUT:
        json.dump([data, bulge_data], OUT, indent=2)
