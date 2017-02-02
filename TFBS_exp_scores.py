import numpy as np
import pybedtools as pbt
import sys
from tqdm import tqdm

def bed_to_chrom_wise_interval(bed, chroms):
    chrom_wise_peaks = {x: [] for x in chroms}
    for i in bed:
        chrom_wise_peaks[i.chrom].append((i.start, i.end))
    return chrom_wise_peaks

if __name__ == "__main__":
    bed_file = sys.argv[1]
    out_file_prefix = sys.argv[2]

    offset=50
    exp_dir='/home/parashar/scratch/quadcomb/data/QG_expectation_scores/'
    chroms=['chr'+str(x) for x in range(1,23)] + ['chrX', 'chrY']
    chrom_info_file='/home/parashar/scratch/hg19_resource/hg19.genome'

    chrom_wise_peaks = bed_to_chrom_wise_interval(pbt.BedTool(bed_file), chroms)
    interval_exp = []
    shuffle_exp = []
    for chrom in chroms:
        pos = np.load("%s/%s_positive.npy" % (exp_dir, chrom), mmap_mode='r')
        neg = np.load("%s/%s_negative.npy" % (exp_dir, chrom), mmap_mode='r')
        
        for i in chrom_wise_peaks[chrom]:
            peak_mid = i[0] + (i[1] - i[0])
            interval_exp.append(pos[peak_mid - offset : peak_mid + offset] +
                                neg[peak_mid - offset : peak_mid + offset])
    
        for i in tqdm(range(100), desc=chrom):
            chrom_wise_peaks = bed_to_chrom_wise_interval(
                pbt.BedTool(bed_file).shuffle(chrom=True, g=chrom_info_file), chroms)
            temp = []
            for i in chrom_wise_peaks[chrom]:
                peak_mid = i[0] + (i[1] - i[0])
                temp.append(pos[peak_mid - offset : peak_mid + offset] +
                            neg[peak_mid - offset : peak_mid + offset])
            shuffle_exp.append(temp)

    np.save("%s_peaks_exp_scores" % out_file_prefix, np.array(interval_exp))
    np.save("%s_shuffle_exp_scores" % out_file_prefix, np.array(shuffle_exp))
