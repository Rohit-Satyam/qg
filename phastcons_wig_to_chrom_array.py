import sys
import numpy as np

def make_array(c):
    chromInfo = {
        '1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276,
        '5': 180915260, '6': 171115067, '7': 159138663, 'X': 155270560,
        '8': 146364022, '9': 141213431, '10': 135534747, '11': 135006516,
        '12': 133851895, '13': 115169878, '14': 107349540, '15': 102531392,
        '16': 90354753, '17': 81195210, '18': 78077248, '20': 63025520,
        'Y': 59373566, '19': 59128983, '22': 51304566, '21': 48129895
    }
    return np.zeros(chromInfo[c])

def read_wiggle(wig):
    lctr = 0
    n = 0
    with open(wig) as h:
        for l in h:
            d = l.rstrip('\n').split(' ')
            if len(d) > 1:
                start_pos = int(d[2].split('=')[1]) - 1
                lctr+=n
                print ('\r%d' % lctr, end='')
                n = 0
            else:
                yield (start_pos+n, float(d[0]))
                n+=1

if __name__ == "__main__":
    wig_dir = '/home/parashar/scratch/quadcomb/data/phastcons/downloaded_wiggle'
    out_dir = '/home/parashar/scratch/quadcomb/data/phastcons/chrom_arrays'
    chrom = sys.argv[1]

    chrom_array = make_array(chrom)
    wig_gen = read_wiggle("%s/chr%s.phastCons46way.primates.wigFix" % (wig_dir, chrom))
    for i in wig_gen:
        chrom_array[i[0]] = i[1]

    print ('Saving array')
    np.save("%s/chr%s" % (out_dir, chrom), chrom_array)