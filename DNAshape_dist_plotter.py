import h5py
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter
import sys


def load_array(fn, name):
    f = h5py.File(fn, 'r')
    a = np.nan_to_num(np.array(f.get(name), dtype=float).flatten())
    f.close()
    return a


def get_value_counts(counter):
    k, v = [], []
    for i in val_count.items():
        k.append(i[0])
        v.append(i[1])
    sorted_idx = sorted(range(len(k)), key=lambda x: k[x])  # sort by keys
    k = np.array(k)[sorted_idx]
    v = np.array(v)[sorted_idx]
    return k, v


def label(ax, shape_type, chrom):
    if shape_type == 'MGW':
        ax.set_xlabel('Minor Groove Width (MGW) in Å')
        ax.set_ylabel('Frequency of positions\nin the chromosome (x1K)')
        ax.set_title('Distribution of minor groove width in %s' % chrom)
    elif shape_type == 'Roll':
        ax.set_xlabel('Roll °')
        ax.set_ylabel('Frequency of positions\nin the chromosome (x1K)')
        ax.set_title('Distribution of degree of DNA roll in %s' % chrom)
    elif shape_type == 'HelT':
        ax.set_xlabel('Helix twist (HelT) °')
        ax.set_ylabel('Frequency of positions\nin the chromosome (x1K)')
        ax.set_title('Distribution of degree of helix twist in %s' % chrom)
    elif shape_type == 'ProT':
        ax.set_xlabel('Propeller twist (ProT) °')
        ax.set_ylabel('Frequency of positions\nin the chromosome (x1K)')
        ax.set_title('Distribution of degree of propeller twist in %s' % chrom)
    else:
        return False


def clean_axis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


if __name__ == '__main__':
    hdf5_file = sys.argv[1]
    shape_type = hdf5_file.split('/')[-1].split('.')[0].split('_')[1]
    chrom = hdf5_file.split('/')[-1].split('.')[0].split('_')[0]

    array = load_array(hdf5_file, 'shape/%s' % shape_type)

    val_count = dict(Counter(array))  # cont frequency of each value
    val_count.pop(0, None)  # remove 'zero' key they will be pain when plotting

    vals, frequency = get_value_counts(val_count)
    v_cumsum = np.cumsum(vals)

    plt.style.use('seaborn-poster')
    fig, ax = plt.subplots(1, 1, figsize=(12, 4))
    _ = ax.bar(vals, frequency / 1000, width=0.02,
               edgecolor='none', color='dodgerblue')
    ax.axvline(np.median(list(val_count.keys())),
               ls='--', c='crimson', lw=2, label='Mean')
    ax.text(0.05, 0.9, 'Total positions=%d' % np.sum(list(val_count.values())),
            ha='left', va='top', transform=ax.transAxes, fontsize=16)
    label(ax, shape_type, chrom)
    ax.legend()
    clean_axis(ax)
    fig.tight_layout()
    plt.savefig('%s_%s.tiff' % (chrom, shape_type), dpi=150)
    plt.close()
