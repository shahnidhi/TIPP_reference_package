*.fna - nucleotide sequences for the 40 marker genes
*.faa - amino acid sequences for the 40 marker genes
Note that the sequence id in the *.fna and the corresponding .faa file match. 
report.txt has the summary information for these sequences combined
	1. sequence id
	2. accession from which sequence was extracted
	3. protein id (extracted from refseq annotation)
	4. length of the protein sequence (aa)
	5. length of the nucleotide sequence (bp)
	6. flag that there is a match in the nucleotide and the corresponding amino acid sequence (1 if nuc_seq_len == aa_seq_len*3 + 3)
	7. taxid of the accession (organism)
	8. ftp path for the accession
	9. original header for the sequence in the cds file from refseq ftp
	10. flag that the accession is from bacteria domain or archaea (bacteria = 1, archaea = 0) 

genes_per_acc.txt is a tsv file with columns
	1. accession
	2. num_genes
	3. refseq_taxid 
	4. species_taxid
	5. organism_name
	6-12: Taxonomy (Kingdom,Phylum,Class,Order,Family,Genus,Species)
	13:  assembly_level (This is just so that we can check accessions that were complete genome should have ~40 marker genes predicted); fragmented/incomplete genomes can have < 40 genes predicted. 

