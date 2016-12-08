
import os
from chimera import runCommand as rc

pdbs = {'2m27': 'VEGF', '4fxm': 'Telomeric', '4wo3': 'cKIT', '5i2v': 'KRAS'}
for i in pdbs:
    os.system('rm ./images/%s.png' % pdbs[i])
    rc("open ./pdb_files/%s.pdb" % i)
    rc("background solid #eeeeee") # matching background color to notebook background
    rc("focus")
    rc("copy file ./images/%s_hd.png width 600 height 600" % pdbs[i])
    rc("close all")