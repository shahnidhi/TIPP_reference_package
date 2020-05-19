import os
import sys

# open RplB_COG0090.species.old_new.mapping
old2new = {}
with open(sys.argv[1]) as f:
	for line in f:
		if line.startswith('tax_id'):
			continue
		else:
			val = line.strip().split(',')
			old2new[val[0]] = val[1]

fw = open(sys.argv[3], 'w')
fw.write('seqname,tax_id\n')
# open RplB_COG0090.species.mapping 
with open(sys.argv[2]) as f:
	for line in f:
		val = line.strip().split(',')
		if val[1] in old2new:
			fw.write(val[0]+','+old2new[val[1]]+'\n')

fw.close()
# Usage python update_species_mapping.py species.old2new.mapping species.mapping species.updated.mapping
# This step should be done only if the taxids have changed in the latest NCBI database
