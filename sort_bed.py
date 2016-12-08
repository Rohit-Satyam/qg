import sys
import os

bed_file = sys.argv[1]

os.system("sort -k 1,1 -k 2,2n %s > %s"  % (bed_file, bed_file+'_temp'))
os.system('rm %s' % bed_file)
os.system('mv %s %s' % (bed_file+'_temp', bed_file))
