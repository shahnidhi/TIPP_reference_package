import os
import sys
from glob import glob

from itertools import groupby

def fasta_iter(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip()
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq
    fh.close()

dirname = sys.argv[1]
files = glob(dirname+'/*faa')
fw = open(sys.argv[2], 'w')
for f in files:
	name = f.split('/')[-1].split('.')[0]
	fiter = fasta_iter(f)
	for ff in fiter:
		fw.write(ff[0] +'\t'+name+'\n')
	

