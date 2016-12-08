import os

samples = {
    'Na_K_1': ['SRR1693705', 'SRR1693706'],
    'Na_K_2': ['SRR1693707', 'SRR1693708'],
    'Na_PDS_1': ['SRR1693709', 'SRR1693710'],
    'Na_PDS_2': ['SRR1693711', 'SRR1693712'],
}

for sample in samples:
    print (sample)
    sub_dirs = samples[sample]
    validated_dir = '../data/g4_seq/reanalysis/oq_validated_bed/%s' % (
        sample)
    d1, d2 = ["/".join([validated_dir, x]) for x in sub_dirs]

    out_switchpoint_dir = '../data/oq_switchpoints'
    oq_bed = out_switchpoint_dir + '/' + sample + '.bed'
    oq_bed_unique = out_switchpoint_dir + '/' + sample + '_unique.bed'
    oq_bed_50 = out_switchpoint_dir + '/' + sample + '_slopped_50.bed'
    oq_fasta = out_switchpoint_dir + '/' + sample + '_slopped_50.fasta'

    os.system("cat %s/*.bed %s/*.bed | sort -k1,1 -k2,2n > %s" % (
        d1, d2, oq_bed))

    OUT = open(oq_bed_unique, 'w')
    with open(oq_bed) as h:
        line = next(h)
        c = line.rstrip('\n').split('\t')
        prev = c[0], c[1], c[5]
        OUT.write(line)
        for line in h:
            c = line.rstrip('\n').split('\t')
            if c[1] == prev[1] and c[5] == prev[2] and c[0] == prev[0]:
                pass
            else:
                OUT.write(line)
            prev = c[0], c[1], c[5]
    OUT.close()

    chrom_info = '~/scratch/hg19_resource/hg19.genome'
    os.system("bedtools slop -i %s -g %s -b 50 > %s" %
              (oq_bed_unique, chrom_info, oq_bed_50))
    genome_file = '~/scratch/hg19_resource/genome.fa'
    os.system("bedtools getfasta -s -fi %s -bed %s > %s" % (
        genome_file, oq_bed_50, oq_fasta))

    temp_intersect_bed = '../data/QuadGraph_%s_intersect_temp.bed' % sample
    os.system("bedtools intersect -S -wao -a %s -b %s > %s" % (
        "../data/QuadGraph.bed", oq_bed, temp_intersect_bed))
    intersect_bed = '../data/oq_switchpoints/qg_oq_%s_intersect.bed' % sample
    OUT = open(intersect_bed, 'w')
    with open(temp_intersect_bed) as h:
        for l in h:
            c = l.rstrip('\n').split('\t')
            if c[-1] == '1':
                OUT.write("\t".join(
                    [c[0], c[1], c[2], c[3], str(int(c[7]) - int(c[1])), c[5]]
                ) + '\n')
    OUT.close()
    os.system("rm %s" % temp_intersect_bed)
    os.system("mkdir ../data/oq_switchpoints/%s" % sample)
    os.system("split -l 500 %s ../data/oq_switchpoints/%s/%s_" % (
        intersect_bed, sample, sample))
