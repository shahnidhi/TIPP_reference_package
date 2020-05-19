import os
import sys
from itertools import groupby

def fasta_iter(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip()
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq
    fh.close()

taxa = {}

with open('genes_per_acc.txt') as f:
	for line in f:
		if line.startswith('#'):
			continue
		val = line.strip().split('\t')
		taxa[val[0]] = val[2]

fw = open(sys.argv[2], 'w')
fw.write('seqname,tax_id\n')
fiter = fasta_iter(sys.argv[1])
for ff in fiter:
	name = '_'.join(ff[0].split('_')[0:2])
	if name in taxa:
		fw.write(ff[0].strip()+','+str(taxa[name])+'\n')
	else:
		print (sys.argv[1]+'\t'+ ff[0]+'\t'+name)

#Usage python generate_species_mapping.py gene_nuc_alignment.fna species.mapping

#cut -f 2 -d ',' RplB_COG0090.species.mapping > RplB_COG0090.species.txt
#taxit update_taxids RplB_COG0090.species.txt -o RplB_COG0090.species.updated.txt
#paste -d "," RplB_COG0090.species.txt <( sed 's/"//g' RplB_COG0090.species.updated.txt) > RplB_COG0090.species.old_new.mapping
#python update_species_mapping.py  RplB_COG0090.species.old_new.mapping RplB_COG0090.species.mapping RplB_COG0090.updated.species.mapping
#taxit taxtable -i RplB_COG0090.species.updated.txt -o RplB_COG0090.taxonomy.table
#python build_taxonomic_tree.py RplB_COG0090.taxonomy.table RplB_COG0090.species.updated.txt RplB_COG0090.unrefined.taxonomy
#perl build_unrefined_tree.pl RplB_COG0090.updated.species.mapping RplB_COG0090.unrefined.taxonomy RplB_COG0090.unrefined.taxonomy.renamed

