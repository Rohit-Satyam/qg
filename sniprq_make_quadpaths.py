import sys
import h5py
from collections import deque
from SNIPRQ import StemGroups, QuadPaths
import numpy as np

if __name__ == "__main__":

    chrom = sys.argv[1]
    base = sys.argv[2]
    quadpaths_fn = sys.argv[3]
    strand = 'positive' if base == 'G' else 'negative'
    max_looplen = 25
    max_bulge_len = 5
    
    out_fn = h5py.File("%s.hdf5" % quadpaths_fn, "w")
    fasta_loc = "/home/parashar/scratch/hg19_resource/chromosomes"
    chrom_seq = "".join([x.rstrip('\n').upper() for x in 
                         open('%s/%s.fa' % (fasta_loc, chrom)).readlines()[1:]])
    seq_qps = deque()
    seq_infos = {
        'graph_end': deque(),
        'graph_start': deque(),
        'num_edges': deque(),
        'num_nodes': deque(),
        'num_quadpaths': deque()
    }
    sg_obj = StemGroups(sequence=chrom_seq, max_loop_len=max_looplen,
                        max_bulge_len=max_bulge_len, base=base)
    for sg in sg_obj.generate_stems():
        qp_obj = QuadPaths(stems=sg, max_loop_len=max_looplen)
        if qp_obj.QPS.shape[0] > 0:
            seq_qps.append(qp_obj.QPS)
            out_fn.create_dataset(str(qp_obj.QPS[0][0][0]), data= qp_obj.QPS)
            for k in qp_obj.info:
                seq_infos[k].append(qp_obj.info[k])
    out_fn.close()
    
    seq_infos = np.array([seq_infos['num_nodes'], seq_infos['num_edges'],
                 seq_infos['num_quadpaths'], seq_infos['graph_start'], seq_infos['graph_end']])
    np.save("%s.npy" % quadpaths_fn, seq_infos)
