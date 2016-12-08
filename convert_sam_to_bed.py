import sys
import os

# Script designed for G4seq.
# Few peculariteis to bear in mind are:
# - Saves mapping quality in the numeric value column of the BED file.
# - Expects a '.' in the reference name and splits it to get the second element

if __name__ == "__main__":
    sam = sys.argv[1]
    out_dir = sys.argv[2].rstrip('/')

    bed = []
    with open(sam) as h:
        while True:
            if next(h)[0:3] == '@PG':
                break
        for l in h:
            c = l.split('\t')
            if c[2] == "*":
                continue
            ref = c[0].split('.')[1]  # <-
            strand = '+' if c[1] == '0' else '-'
            bed.append(
                "\t".join(
                    [c[2], str(int(c[3]) - 1), c[3], ref, c[4], strand]
                )
            )

    out_name = "%s/%s.bed" % (out_dir, sam.split('/')[-1].split('.')[0])
    with open(out_name, 'w') as OUT:
        OUT.write("\n".join(bed))
    os.system("mv %s %s" % (out_name, out_name + '_temp'))
    os.system("sort -k 1,1 -k 2,2n %s > %s" % (out_name + '_temp', out_name))
    os.system("rm %s" % out_name + '_temp')
