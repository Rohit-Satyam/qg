import numpy as np
import glob
from scipy.stats import poisson, norm
import json

base_dir = "../data/ENCODE_TFBS_Uniform/exp_scores"
files = glob.glob("%s/*peaks_exp_scores.npy" % base_dir)

samples = {}
for fn in files:
    try:
        cell_line, antibody = fn.split('/')[-1].split('_peaks_exp_scores.npy')[0].split('_')
    except ValueError:
        cell_line, _, antibody = fn.split('/')[-1].split('_peaks_exp_scores.npy')[0].split('_')
    if antibody not in samples:
        samples[antibody] = {}
    samples[antibody][cell_line] = fn
    
def get_pvals(peak_array, shuffle_array):
    reshaped_shuffle_array = [[] for x in range(100)]
    for i in range(0, 2400, 100):
        for n,j in enumerate(shuffle_array[i:i+100]):
            reshaped_shuffle_array[n].extend(j)
    reshaped_shuffle_array = np.array(reshaped_shuffle_array)
    if len(reshaped_shuffle_array.shape) == 3:
        try:
            shuffle_sum = [reshaped_shuffle_array[x].sum() for x in range(100)]
        except ValueError:
            return {}
        poisson_pval = poisson.pmf(peak_array.sum(), mu=np.mean(shuffle_sum))*len(peak_array)
        norm_pval = norm.pdf(peak_array.sum(), loc=np.mean(shuffle_sum),
                         scale=np.std(shuffle_sum))*len(peak_array)
        return {
            'p_pval': poisson_pval,
            'n_pval': norm_pval,
            'tot_peaks': len(peak_array),
            'sum_peaks': peak_array.sum(),
            'sum_shuffle': np.mean(shuffle_sum),
            'std_shuffle': np.std(shuffle_sum)
        }
    return {}

results = {}
for antibody in samples:
    results[antibody] = {}
    for cell_line in samples[antibody]:
        print (cell_line, antibody)
        fn = samples[antibody][cell_line]
        peak_array = np.load(fn)
        shuffle_array = np.load(fn.replace('peaks', 'shuffle'))
        results[antibody][cell_line] = get_pvals(peak_array, shuffle_array)

with open("../data/ENCODE_TFBS_Uniform/pvals.json", 'w') as OUT:
    json.dump(results, OUT, indent=2)