import json
from tqdm import tqdm
import os

def get_skip_lines(fn):
    with open(fn) as h:
        for n,l in enumerate(h):
            if l[0] != '#':
                break
    return n

def gen_transcript_line(fn, skip_num):
    with open(fn) as h:
        for i in range(skip_num):
            next(h)
        for l in h:
            c = l.rstrip('\n').split('\t')
            if c[2] == 'transcript':
                yield c

def get_transcript_id_gene_name(split_line):
    i = split_line[8].split(';')
    tid = i[1].strip(' ').split(' ')
    if tid[0] == 'transcript_id':
        tid = tid[1].strip('"')
    else:
        print ('WARNING: "transcript_id" not found in position 1')
    g_name = i[4].strip(' ').split(' ')
    if g_name[0] == 'gene_name':
        g_name = g_name[1].strip('"')
    else:
        print ('\nWARNING: "gene_name" not found in position 4')
    return (tid, g_name)

def get_tss(split_line):
    if split_line[6] == '+':
        tss = int(split_line[3]) - 1 # GTF format is 1 based
    elif split_line[6] == '-':
        tss = int(split_line[4]) - 1 # End is inclusive. otherwise would have been -2
    else:
        print ('\nWARNING: Unknown strand symbol detected')
    return tss

if __name__ == "__main__":
    gtf_fn = '/home/parashar/scratch/quadcomb/data/annotation/gencode.v25lift37.annotation.gtf'
    out_bed = '/home/parashar/scratch/quadcomb/data/annotation/gencode_tss.bed'
    out_map = '/home/parashar/scratch/quadcomb/data/annotation/gencode_tid_to_gene_name.json'
    file_num_lines = 200139

    num_line_skip = get_skip_lines(gtf_fn)
    bed = []
    t_id_gene_map = {}
    line_generator = gen_transcript_line(gtf_fn, num_line_skip)

    for i in tqdm(line_generator, total=file_num_lines):
        tid, gene = get_transcript_id_gene_name(i)
        tss = get_tss(i)
        if tid not in t_id_gene_map:
            t_id_gene_map[tid] = gene
        else:
            print ('\nWARNING: entry attempted with duplicate "transcript_id"')
        bed.append('\t'.join([i[0], str(tss), str(tss+1), tid, '.', i[6]]))
    print ('\n')
    with open(out_map, 'w') as OUT:
        json.dump(t_id_gene_map, OUT, indent=2)
    with open('temp_gencode_tss.bed', 'w') as OUT:
        OUT.write('\n'.join(bed))
    print ("Sorting BED file")
    os.system("bedtools sort -chrThenSizeA -i temp_gencode_tss.bed > %s" % out_bed)
    os.system('rm temp_gencode_tss.bed')
