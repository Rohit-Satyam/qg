import sys
import os

fastq_file = sys.argv[1]
chunk_lines = int(sys.argv[2])
out_dir = sys.argv[3]

if chunk_lines % 4 != 0:
    raise ValueError('Chunk lines should be multiple of 4')

file_ctr = 0

dir1 = os.path.join(out_dir, 'control')
dir2 = os.path.join(out_dir, 'treated')
if os.path.isdir(dir1):
    os.system("rm -rf %s" % dir1)
if os.path.isdir(dir2):
    os.system("rm -rf %s" % dir2)
os.system("mkdir -p %s" % dir1)
os.system("mkdir -p %s" % dir2)

OUT1 = open(os.path.join(dir1, '%d.fastq' % file_ctr), 'w')
OUT2 = open(os.path.join(dir2, '%d.fastq' % file_ctr), 'w')

with open(fastq_file) as handle:
    for n, line in enumerate(handle):
        OUT1.write(line)
        OUT2.write(line)
        line = next(handle)
        OUT1.write(line[:151] + '\n')
        OUT2.write(line[151:])
        line = next(handle)
        OUT1.write(line)
        OUT2.write(line)
        line = next(handle)
        OUT1.write(line[:151] + '\n')
        OUT2.write(line[151:])
        x = n + 1
        if x >= chunk_lines and x % chunk_lines == 0:
            OUT1.close()
            OUT2.close()
            file_ctr += 1
            OUT1 = open(os.path.join(out_dir, 'control', '%d.fastq' % file_ctr), 'w')
            OUT2 = open(os.path.join(out_dir, 'treated', '%d.fastq' % file_ctr), 'w')
OUT1.close()
OUT2.close()
