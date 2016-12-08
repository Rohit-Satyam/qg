import glob
import json
from tqdm import tqdm
import sys


chrom = sys.argv[1]
strand = sys.argv[2]
out_dir = '/home/parashar/scratch/quadcomb/data/chrom_wise_graph_summary/%s' % strand

jsons = glob.glob('/home/parashar/scratch/quadcomb/data/quad_graphs_mll50/%s/%s/*/*.json' % (
                  chrom, strand))

data = {
    'nodes': [],
    'edges': [],
    'quads': [],
    'beginPos': [],
    'ghostBeginPos': [],  # begin pos before nodes were deleted
    'endPos': []
}
for j in tqdm(jsons):
    d = json.load(open(j))
    for i in d:
        data[i].append(d[i])
    data['ghostBeginPos'].append(int(j.split('/')[-1].split('.')[0]))

with open('%s/%s.json' % (out_dir, chrom), 'w') as OUT:
    json.dump(data, OUT, indent=2)
