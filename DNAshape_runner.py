import os

base_dir = '/home/parashar/scratch/quadcomb/data/dna_shape'
chrom_dir = '/home/parashar/scratch/hg19_resource/chromosomes'
save_dir = 'hdf5_files'
script_dir = 'dnashapeR_scripts'

base_file = """
library(DNAshapeR)
library(rhdf5)

pred <- getShape('%s/%s.fa', shapeType='All', parse=TRUE)
h5createFile("%s_MGW.h5")
h5createFile("%s_HelT.h5")
h5createFile("%s_ProT.h5")
h5createFile("%s_Roll.h5")
h5createGroup("%s_MGW.h5", "shape")
h5createGroup("%s_HelT.h5", "shape")
h5createGroup("%s_ProT.h5", "shape")
h5createGroup("%s_Roll.h5", "shape")
h5write(pred$MGW, file="%s_MGW.h5", name="shape/MGW")
h5write(pred$HelT, file="%s_HelT.h5", name="shape/HelT")
h5write(pred$ProT, file="%s_ProT.h5", name="shape/ProT")
h5write(pred$Roll, file="%s_Roll.h5", name="shape/Roll")
"""

# if os.path.isdir('%s/%s' % (base_dir, save_dir)):
#     os.system('rm -rf %s/%s' % (base_dir, save_dir))
# if os.path.isdir('%s/%s' % (base_dir, script_dir)):
#     os.system('rm -rf %s/%s' % (base_dir, script_dir))
# os.system('mkdir %s/%s' % (base_dir, save_dir))
# os.system('mkdir %s/%s' % (base_dir, script_dir))

chroms = ['chr' + str(x) for x in range(15, 23)] + ['chrX', 'chrY']

for chrom in chroms:
    fn_path = "%s/%s/%s" % (base_dir, save_dir, chrom)
    arg_list = [fn_path for i in range(12)]
    file_content = base_file % (chrom_dir, chrom, *arg_list)
    script_name = '%s/%s/%s.Rscript' % (base_dir, script_dir, chrom)
    with open(script_name, 'w') as h:
        h.write(file_content)
    #cmd = "bsub -q debugq -n 4 -J %s -R 'span[hosts=1]' -o %s.log Rscript %s" % (chrom, chrom, script_name)
    cmd = "Rscript %s" % script_name
    os.system(cmd)
