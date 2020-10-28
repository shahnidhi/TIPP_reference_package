# TIPP reference packages 

## Quick Links
- TIPP paper - [https://doi.org/10.1093/bioinformatics/btu721](https://doi.org/10.1093/bioinformatics/btu721)
- TIPP software can be found [here](https://github.com/smirarab/sepp/blob/master/README.TIPP.md).
- TIPP reference dataset (2014) - https://github.com/tandyw/tipp-reference/releases/download/v2.0.0/tipp.zip
- **(New)** TIPP reference dataset (2020) - https://obj.umiacs.umd.edu/tipp/tipp2-refpkg.tar.gz 
<!---
- **(New)** Sparse TIPP reference dataset 1 (2020); contains only one sequence per species - ADD LINK
- **(New)** Sparse TIPP reference dataset 1 (2020); contains only two sequences per genus - ADD LINK
--->

In this document, we describe the protocol used to construct a new version of [TIPP](https://doi.org/10.1093/bioinformatics/btu721) reference packages. We used the same set of 40 marker genes as used by Mende et al. 2013, Sunagawa et al. 2013, and mOTUs. These marker genes are believed to be single-copy and universally present in prokaryotic genomes. 

## Data sources
We downloaded all Bacterial and Archael genomes from the NCBI RefSeq database. RefSeq provides a metadata file for both Bacteria and Archaea genomes. This file contains useful information such as genome accession, taxid, name, ftp download information, etc. We downloaded these files in November 2019,, and are provided in the data folder. One can download the latest version of files from
- `wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt`
- `wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt`

## Extracting gene sequences and group them by markers
For each genome accession in this master list, we download the genome sequence data, protein sequences, and nucleotide gene sequences. We modify the names of the protein and nucleotide gene sequences such that the corresponding protein and nucleotide gene sequence have same name. We also add genome accession to the start of the gene sequence, to make it possible to track the origin of the gene sequence in the database.
For example, one of the gene sequence in ArgS (COG0018) marker gene is `GCF_002287175.1_NZ_LMVM01000012.1_cds_WP_069582217.1_1177`, which has genome accession, sequence accession, and the protein name information in it's identifier. 
Once, we have protein and corresponding nucleotide gene sequences from a genome, we used [fetchMG](http://vm-lux.embl.de/~mende/fetchMG/about.html) tool to extract 40 marker gene sequences. Note that fetchMG uses both protein and nucleotide sequences, so keep the name of the two files same, just change the extension to .faa for proteins ana .fna for nucleotide sequences. We used fetchMG v1.0 for our analysis. \
```bash
fetchMG.pl -m extraction -v gene_aa_filename -o output_folder
```

The script get_sequences.py combines the steps of downloading genomes and extracting 40 marker genes, just run
```bash
python get_sequences.py
```
This should create an output folder with \*.fna (nucleotide gene sequences) and \*.faa (protein sequences) files for the 40 marker genes. The report.txt has metadata information for the selected gene sequences. Please look at the [output_README.txt](https://github.com/shahnidhi/TIPP_reference_package/blob/master/output_README.txt) file for column headers.  

## Filtering gene sequences
We removed all gene sequences that had a mismatch in their nucleotide and protein sequence lengths. The nucleotide gene sequence length (after removing stop codon) should match 3\*protein sequence length. We also removed sequences that were 3 std deviation away from the median gene length. There are ~140K-170K gene sequences per marker gene after filtering based on length. 

## Alignment Estimation
For each marker gene, we generate multiple sequence alignment (MSA) of the protein sequences using UPP software. The parameters chosen essentially generates PASTA alignments, because we have carefully chosen gene sequences such that they are almost all full-length. We used [UPP](https://github.com/smirarab/sepp/blob/master/README.UPP.md) version 4.3.10 software.

```bash
run_upp.py -s gene_aa_filename.faa -p tmp_dir  -B 1000000 -M -1 \
-T 0.33 -m amino -o output_folder
```

We used alignment with insertion sites masked (\*masked.fasta) in all subsequent steps.
### Translate back to DNA
Once we generate protein MSAs, we translate the proteins back to nucleotide sequences. Because we have nucleotide sequences for each protein sequence, we use that information to match amino acid to the corresponding codon. 
```bash
python backtranslate_refseq.py protein_alignment.fasta \
       gene_nuc_filename.fna gene_nuc_alignment.fasta output_log.txt
```
### Remove gappy sites and fragmentary data
We removed gappy sites from the alignment. This step is performed after translating back to DNA, so that we don't lose amino acid to codon mapping from the RefSeq files. 
We removed all sites with 95% or more gaps.
```bash
f=[alignment file]
percent=5
m=`echo $( grep ">" $f|wc -l ) \* $percent / 100 |bc`
run_seqtools.py -infile $f -masksites $m -outfile $f.mask${percent}sites.fasta
```
We also remove sequences that are fragments.
```bash
taxapercent=33 
m2=`echo $( cat $f.mask${percent}sites.fasta|wc -L ) \* $taxapercent / 100 |bc`
run_seqtools.py -infile $f.mask${percent}sites.fasta -filterfragments $m2 -outfile $out
```
## Gene Tree Estimation and Refining The Taxonomy
### Generate Unrefined taxonomy
Based on Refseq metadata file, we have taxid for each genome sequence. However, sometimes NCBI taxid gets depreciated, merged, or updated. To get the latest taxid and the complete NCBI lineage for each gene sequence, we run the following steps. We rely heavily on [Taxtastic](http://fhcrc.github.io/taxtastic/) v0.8.11 software. 
```bash
python generate_species_mapping.py gene_nuc_alignment.fna species.mapping
cut -f 2 -d ',' species.mapping > species.txt
taxit update_taxids species.txt -o species.updated.txt
taxit taxtable -i species.updated.txt -o taxonomy.table

paste -d "," species.txt <( sed 's/"//g' species.updated.txt) > species.old2new.mapping
python update_species_mapping.py species.old2new.mapping \ 
        species.mapping species.updated.mapping
        
python build_taxonomic_tree.py taxonomy.table species.updated.txt unrefined.taxonomy
perl build_unrefined_tree.pl species.updated.mapping \ 
     unrefined.taxonomy unrefined.taxonomy.renamed
```
The output of “taxit taxtable” is a table with all the taxonomic ranks organized. This will be later used in computing taxonomic profile. 
### Generate Refined taxonomy
We used GTRCAT and GTRGAMMA model of [RAxML](https://github.com/stamatak/standard-RAxML) (version 8.2.12) to generate refined taxonomy for each marker gene. 
```bash
raxmlHPC-PTHREADS-AVX -j -m GTRCAT -F -T 4 -p 1111 -g unrefined.taxonomy.renamed \
           -s gene_nuc_alignment.fasta -n refined -w ${work}/raxml_output/

raxmlHPC-PTHREADS-AVX -j -m GTRGAMMA -f e -t RAxML_result.refined -T 4 -p 1111 \
           -s gene_nuc_alignment.fasta -n optimized -w ${work}/raxml_output
```
### Generate ML gene trees
We also generate a maximum likelihood gene tree for each marker gene. 
```bash
raxmlHPC-PTHREADS-AVX -j -m GTRCAT -F -T 4 -p 1111 \ 
           -s gene_nuc_alignment.fasta -n mlgene -w {work}/raxml_output_mlgene/

raxmlHPC-PTHREADS-AVX -j -m GTRGAMMA -f e -t RAxML_result.mlgene -T 4 -p 1111 \ 
         -s gene_nuc_alignment.fasta -n optimized -w ${work}/raxml_output_mlgene/
```
## Databases: BLAST
We concatenate all gene sequences (nucleotide) from all marker genes to create a combined fasta file, and a sequence to marker mapping file (seq2marker.tab).
We used [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (v 2.9.0) to create database files.
```bash
    cat <all nucleotide gene sequence files (*.fna)> > alignment.fasta
    makeblastdb -in alignment.fasta -out alignment.fasta -dbtype nucl
```


