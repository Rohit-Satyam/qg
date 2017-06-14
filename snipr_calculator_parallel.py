import sys
import numpy as np
import multiprocessing as mp
import time
from SNIPRQ import QuadPaths, StemGroups
from collections import deque

def load_table():
    look_up_table = {}
    with open('../data/table_stem_loop_stem_frequency.csv') as h:
        next(h)
        for l in h:
            c = l.rstrip('\n').split(',')
            look_up_table[(int(c[0]), int(c[1]), int(c[2]))] = int(c[5])
            look_up_table[(int(c[2]), int(c[1]), int(c[0]))] = int(c[5])
    return look_up_table

def scorer_wrapper(sg):
    
    def scorer(qp):
        loop_lens = qp[:, 0][1:] - qp[: , 2][:-1] - 1
        bulge_lens = [x[2] - x[0] - 2 for x in qp]
        a = 1/look_up_table[(bulge_lens[0], loop_lens[0], bulge_lens[1])]
        b = 1/look_up_table[(bulge_lens[1], loop_lens[1], bulge_lens[2])]
        c = 1/look_up_table[(bulge_lens[2], loop_lens[2], bulge_lens[3])]
        return int(3/(a+b+c))

    start = sg[0][0]
    stop = sg[-1][-1]
    max_scores = np.zeros(stop-start+1, dtype=int)
    qp_obj = QuadPaths(stems=sg, max_loop_len=25)
    for qp in qp_obj.QPS:
        score = scorer(qp)
        qpf = qp.flatten()-start
        max_scores[qpf[max_scores[qpf] < score]] = score
    if qp_obj.info != []:
        return (start, stop+1, max_scores, qp_obj.info)
    return False

if __name__ == '__main__':
    
    chrom = sys.argv[1]
    base = sys.argv[2]
    num_processes = int(sys.argv[3])
    
    max_looplen = 25
    max_bulge_len = 5
    look_up_table = load_table()

    mpPool = mp.Pool(num_processes)
    chrom_seq = ''.join([x.rstrip('\n') for x in 
                     open('../../hg19_resource/chromosomes/%s.fa' %
                          chrom).readlines()[1:]])
    sg_obj = StemGroups(sequence=chrom_seq, max_loop_len=max_looplen,
                        max_bulge_len=max_bulge_len, base=base)
    chrom_scores = np.zeros(len(chrom_seq), dtype=int)
    chrom_info = []
    sg_bucket = deque()
    for sg in sg_obj.generate_stems():
        print ('\r%d' % sg_obj.pos, end='')
        sg_bucket.append(sg)
        if len(sg_bucket) == 5000:
            print ('yoo')
            mpPool = mp.Pool(num_processes)
            res_block = mpPool.map(scorer_wrapper, sg_bucket)
            mpPool.close()
            mpPool.join()
            for res in res_block:
                if res is not False:
                    chrom_scores[res[0]:res[1]] = res[2]
                    chrom_info.append(res[3])
                    sg_bucket = []
    if len(sg_bucket) > 0:
        mpPool = mp.Pool(num_processes)
        res_block = mpPool.map(scorer_wrapper, sg_bucket)
        mpPool.close()
        mpPool.join()
        for res in res_block:
            if res is not False:
                chrom_scores[res[0]:res[1]] = res[2]
                chrom_info.append(res[3])
                sg_bucket = []
    strand = 'positive' if base == 'G' else 'negative'
    np.save('../data/snipr/%s_%s_scores' % (chrom, strand), chrom_scores)
    np.save('../data/snipr/%s_%s_info' % (chrom, strand), chrom_info)
