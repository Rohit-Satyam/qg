import numpy as np
import networkx as nx
from tqdm import tqdm
from collections import deque
import json
import os
import sys
import multiprocessing as mp
import time


def prepare_dirs(saveDir, strand, dirInterval, seq_len):
    if not os.path.isdir(saveDir):
        os.mkdir(saveDir)
    d = os.path.join(saveDir, strand)
    if os.path.isdir(d):
        os.system("rm -rf %s" % d)
    os.mkdir(d)
    for i in range(dirInterval, seq_len + dirInterval, dirInterval):
        os.mkdir(os.path.join(saveDir, strand, str(i)))
    return True


def generate_stems(seq, base, bulge_len, maxLoopLen):
    stemLen = 3  # QuadGraph supports stem length of 3 only. DON'T CHANGE!!
    kmerSize = stemLen + bulge_len
    stems = []
    for p in tqdm(range(len(seq) - kmerSize + 1)):
        s = []
        kmer = np.array(list(seq[p: p + kmerSize]))
        if kmer[0] == base:
            if kmer[1] == base:
                s.extend(([[p, p + 1, p + i] for i in
                           np.where(kmer == base)[0][2:]]))
            s.extend([[p, p + i, p + i + 1] for i in
                      range(2, kmerSize - 1) if
                      kmer[i] == base and kmer[i + 1] == base])
        if len(s) > 0:
            if len(stems) > 1 and s[0][0] - stems[-1][2] > maxLoopLen:
                if len(stems) > 4:
                    yield stems  # releasing stems
                stems = list(s)
            else:
                stems.extend(s)
    for j in range(1, len(kmer) - stemLen):  # loop for last kmer
        xp = p + j
        k = kmer[j:]
        if k[0] == base:
            if k[1] == base:
                stems.extend(([[xp, xp + 1, xp + i] for i in
                               np.where(k == base)[0][2:]]))
            stems.extend([[xp, xp + i, xp + i + 1] for i in
                          range(2, len(k) - 1) if
                          k[i] == base and k[i + 1] == base])
    if len(stems) > 4:
        yield stems


def stem_encoder(stem_array, norm):
    diff = np.diff(stem_array)
    if diff[0] != 1:
        return '%d%s**' % (stem_array[0] - norm, '-' * (diff[0] - 1))
    else:
        return '%d*%s*' % (stem_array[0] - norm, '-' * (diff[1] - 1))


def stem_decoder(code, norm):
    bulge = code.count('-')
    first_val = int(code.replace('*', '').replace('-', '')) + norm
    if code[-2] == '*':
        return [first_val, first_val + bulge + 1, first_val + bulge + 2]
    else:
        return [first_val, first_val + 1, first_val + bulge + 2]


def make_graph(stems, bulgeLen, saveDir, dirInterval, maxLoopLen):
    stemLen = 3  # dont change these values.
    offset = (mll + stemLen + bulgeLen) * (stemLen + bulgeLen + 1)
    stems = np.array(stems)
    dirInterval = int(dirInterval)
    G = nx.DiGraph()
    node_name_norm = stems[0][0]
    for i in range(stems.shape[0] - 1):
        s = stems[i + 1: i + offset][:, 0]
        b = (s > stems[i][2] + 1) & \
            (s <= stems[i][2] + maxLoopLen + 1)
        for j in np.where(b == True)[0]:
            n1 = stem_encoder(stems[i], node_name_norm)
            n2 = stem_encoder(stems[i + j + 1], node_name_norm)
            edge_len = stems[i + j + 1][0] - stems[i][2]
            G.add_edge(n1, n2, length=int(edge_len))
    quads = deque()
    num_quads = 0
    for n1 in G.nodes():
        for n2 in G.successors(n1):
            for n3 in G.successors(n2):
                for n4 in G.successors(n3):
                    quads.extend([n1, n2, n3, n4])
                    num_quads += 1
    del_nodes = list(set(G.nodes()).difference(set(quads)))
    for node in del_nodes:
        G.remove_node(node)
    if len(G.edges()) >= 3 and nx.dag_longest_path_length(G) >= 3:
        decoded_stems = []
        for node in G.nodes():
            decoded_stems.append(stem_decoder(node, node_name_norm))
        decoded_stems = np.array(decoded_stems)
        save_file = ((node_name_norm // dirInterval) *
                     dirInterval) + dirInterval
        save_name = os.path.join(saveDir, str(save_file), str(node_name_norm))
        json_dict = {
            'nodes': len(G.nodes()),
            'edges': len(G.edges()),
            'quads': num_quads,
            'beginPos': int(decoded_stems.min()),  # not same as node_name_norm
            'endPos': int(decoded_stems.max())
        }
        with open('%s.json' % save_name, 'w') as OUT:
            json.dump(json_dict, OUT, indent=2)
        nx.write_graphml(G, '%s.graphml' % save_name)
    return True


if __name__ == "__main__":
    fasta = sys.argv[1]
    out_dir = sys.argv[2]
    base = sys.argv[3]
    bulge_len = int(sys.argv[4])
    mll = int(sys.argv[5])
    num_processes = int(sys.argv[6])
    dirSize = int(1e6)

    seq = open(fasta).readlines()[1:]
    seq = ''.join([x.rstrip('\n') for x in seq])

    if base not in ['G', 'C']:
        print ('base should be either "G" or "C"')
        exit(1)
    strand_dir = 'positive' if base == 'G' else 'negative'

    chrom = fasta.split('/')[-1].split('.')[0]
    prepare_dirs(os.path.join(out_dir, chrom), strand_dir, dirSize, len(seq))
    saveDirectory = os.path.join(out_dir, chrom, strand_dir)

    for stems in generate_stems(seq, base, bulge_len, mll):
        while len(mp.active_children()) >= num_processes:
            time.sleep(0.00001)
        p = mp.Process(target=make_graph, args=(
            stems, bulge_len, saveDirectory, dirSize, mll))
        p.start()
    while len(mp.active_children()) > 0:
        time.sleep(1)
