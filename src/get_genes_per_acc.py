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

count = {}
for f in files:
	fiter = fasta_iter(f)
	seen = {}
	acc = {}
	for ff in fiter:
		accname = '_'.join(ff[0].split('_')[0:2])
		seen[ff[0]] =1
		if accname in acc:
			print ("multiple genes of acc in file", accname, f)
		acc[accname] = 1
	

	fiter = fasta_iter(f.replace('faa','fna'))
	for ff in fiter:
		if ff[0] not in seen:
			print ("print not found in faa", f)
	for accname in acc:
		if accname in count:
			count[accname] += 1
		else:
			count[accname] = 1

## read in taxonomy 
taxa = {}
with open ('taxonomy.info') as f:
	for line in f:
		val = line.strip().split(',')
		taxa[val[0]] = line.strip()

refseq = {}
with open('assembly_summary_refseq.txt') as f:
	for line in f:
		if line.startswith('#'):
			continue
		val =line.strip().split('\t')

		tokeep = [val[5], val[6], val[7], val[11]]
		refseq[val[0]] = tokeep


fw = open(sys.argv[2], 'w')
fw.write('#accession\tnum_genes\trefseq_taxid\tspecies_taxid\torganism_name\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tassembly_level\n')
for accname in count:
	taxaval = taxa[accname].split(',')
	taxaval = [str(x) if x != '' else "NA" for x in taxaval ]
	refseqval = refseq[accname]
	arr2print = [accname, str(count[accname]), str(refseqval[0])] + [str(refseqval[1]), str(refseqval[2])]+ taxaval[2:] + [refseqval[3]]
	fw.write('\t'.join(arr2print)+ '\n')

# Usage python get_genes_per_acc.py  output_folder genes_per_acc.txt
