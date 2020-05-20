import shutil, os
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



def main():
	os.system("mkdir -p output")
	final_fr = open('output/report.txt', 'w')
	final_reprocess = open('output/reprocess_accessions.txt', 'w')
	final_fr.write('#seq_id\taccession\tprotein_name\tlen_gene_nuc\tlen_gene_aa\tnuc_len_match_aa_len\tncbi_taxid\tftp_link\tgene_header_in_refseq\tis_bacteria?\n')
	cognames = {}
	with open('data/FetchMG_COGid2gene_name.tsv') as f:
		for num, line in enumerate(f):
			if num == 0:
				continue
			val = line.strip().split('\t')
			cognames[val[0]] = val[1]
	filehandle_cogs_fna = {}
	filehandle_cogs_faa = {}
	for gene in cognames:
		filehandle_cogs_fna[gene] = open('output/'+cognames[gene]+'_'+gene+'.fna', 'w')
		filehandle_cogs_faa[gene] = open('output/'+cognames[gene]+'_'+gene+'.faa', 'w')

			
	file_accession=['data/assembly_summary_archaea.txt', 'data/assembly_summary_bacteria.txt']
	for bac_arc_ind, refseq_file in enumerate(file_accession):
		with open(refseq_file) as f:
			for num, line in enumerate(f):
				flag2reprocess = False
				val = line.strip().split('\t')
				if line.startswith('#'):
					continue
				accession = val[0]
				taxid = val[5]
				speciestaxid = val[6]
				ftp = val[19]
				trailname = ftp.split('/')[-1]
				#Download protien sequences (aa), and cds sequences (nucleotide) for each genome accession
				cdsfilename = trailname + "_cds_from_genomic.fna"
				proteinfilename = trailname + "_protein.faa"

				cdsftp = ftp + '/' + cdsfilename + ".gz"
				proteinftp = ftp + '/' + proteinfilename + ".gz"
				downloadcds = "wget " + cdsftp	
				downloadprotein = "wget " + proteinftp
				if (os.path.exists(cdsfilename) == False):
					try:
						os.system(downloadcds)
					except Exception as e:
						flag2reprocess = True
						print ("failed downloading cds file for: ", accession) 
					if (os.path.exists(cdsfilename+'.gz') == True):
						try:
							os.system("gunzip -f " + cdsfilename + ".gz")
						except Exception as e:
							print ("failed downloading and unzipping cds file for: ", accession)
							flag2reprocess = True
				if (os.path.exists(proteinfilename) == False):
					try:
						os.system(downloadprotein)
					except Exception as e:
						print ("failed downloading protein file for: ", accession)
						flag2reprocess = True 
					if (os.path.exists(proteinfilename+'.gz') == True):
						try:
							os.system("gunzip -f " + proteinfilename + ".gz")
						except Exception as e:
							print ("failed downloading and unzipping protein file for: ", accession)
							flag2reprocess = True
				if flag2reprocess:
					final_reprocess.write(line.strip()+'\n')
					continue
				proteinseqmap = {}
				fiter = fasta_iter(proteinfilename)
				for ff in fiter:
					header = ff[0].strip().split(' ')[0]
					proteinseqmap[header] = ff[1]

				gene_nuc_filename = accession + "_gene.fna"
				gene_aa_filename = accession + "_gene.faa"
				report_filename = accession + "_gene_summary.txt"
				fgn = open(gene_nuc_filename, 'w')
				fga = open(gene_aa_filename, 'w')
				fgr = open(report_filename, 'w')

				# Keep only those genes that have protein id and are not pseudo genes
				fiternuc = fasta_iter(cdsfilename)
				for ff in fiternuc:
					header = ff[0]
					seq_id = accession + '_' + header.split(' ')[0].split('|')[-1]
					nucseq = ff[1]
					if "protein_id=" in header:
						protein_name = header.split('protein_id=')[1].split(']')[0]
					elif "pseudo=true" not in header:
						print ("#", accession, header, "\tcds with no protein assigned and not pseudo gene")
						continue
					else:
						continue
					if protein_name in proteinseqmap:
						corr_protein = proteinseqmap[protein_name]
					else:
						print ("#", accession, header, "\tprotein missing")
						continue

					nuclen = len(nucseq)
					aalen = len(corr_protein)

					if nuclen - 3  == 3 * aalen:
						matchflag = 1
					else:
						matchflag = 0
						print ("#", accession, header, "\tnucleotide and aa seq length don't match\t", nuclen, aalen, str(nuclen -3 -(aalen*3)), str(bac_arc_ind))
					## Writes nucleotide gene sequences, aa gene sequences and the metadata information for all genes for the accession. Later, we will keep records for only those genes that belong to 40 marker genes. 
					fgr.write('\t'.join([seq_id, accession, protein_name, str(nuclen), str(aalen), str(matchflag), taxid, ftp, header, str(bac_arc_ind)]) + '\n')
					fgn.write('>' + seq_id +'\n' + nucseq +'\n')
					fga.write('>' + seq_id + '\n' + corr_protein + '\n')
				fgr.close()
				fgn.close()
				fga.close()

				## Fetchmg retrieves all genes that align well to one of the 40 marker genes. -v parameter retrieves very best hit for each marker gene. 
				output_folder = accession+"_cogoutput"
				try:
					fetchmg_command = "fetchMG.pl -m extraction -v " + gene_aa_filename + " -o " + output_folder
					os.system(fetchmg_command)
				except Exception as e:
					print ("FetchMG failed", accession)
					continue

				## write to the final files
				genesofinterest = {}
				for gene in cognames:
					## reading nuc sequence first
					if os.path.exists(output_folder+'/'+gene+'.fna'):
						fitergene = fasta_iter(output_folder+'/'+gene+'.fna')
						for ff in fitergene:
							genesofinterest[ff[0]] = 1
							filehandle_cogs_fna[gene].write('>' + ff[0] + '\n' + ff[1] + '\n')

					## reading aa sequence second
					if os.path.exists(output_folder+'/'+gene+'.faa'):
						fitergene = fasta_iter(output_folder+'/'+gene+'.faa')
						for ff in fitergene:
							filehandle_cogs_faa[gene].write('>' + ff[0] + '\n' + ff[1] + '\n')

				## write relevant enteries from report file
				with open(report_filename) as f:
					for line in f:
						val = line.strip().split('\t')
						if val[0] in genesofinterest:
							final_fr.write(line.strip()+'\n')

				##cleanup unnecessary files 
				try:
					os.system("rm " + cdsfilename)
					os.system("rm "+ proteinfilename)
					os.system("rm " + gene_nuc_filename +"*")
					os.system("rm " + gene_aa_filename +"*")
					os.system("rm " + report_filename)
					os.system("cat " + output_folder+"/"+ accession +"_gene.all.marker_genes_scores.table >> " + "output/all.marker_genes_scores.table")
					shutil.rmtree(output_folder)
				except:
					print ("Error while cleaning up", accession)	
					continue
				print ("#completed", accession)

if __name__ == '__main__':
	main()