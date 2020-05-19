# TIPP reference packages 
Here, we describe the procedure used to construct the new version of TIPP reference packages. We used the same set of 40 marker genes as used by Mende et al. 2013, Sunagawa et al. 2013, and mOTUs. These marker genes are believed to be single-copy and universally present in prokaryotic genomes. 

## Data sources
We downloaded all Bacterial and Archael genomes from the NCBI RefSeq database. RefSeq provides a metadata file for both Bacteria and Archaea genomes. This file contains useful information such as genome accession, taxid, name, ftp download information, etc. We downloaded these files in November 2019, and are provided in the data folder. One can download the latest version of files from
- `wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt`
- `wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt`

For each genome accession in this master list, we download the genome sequence data (nuc), protein sequences (aa), and cds sequences (nuc). We modify the names of the protein and cds sequences such that the corresponding protien and cds sequence have same name. We also add genome accession to the start of the gene sequence, to make it possible to track the origin of each gene sequence in the database. 
For example, one of the gene sequence in ArgS (COG0018) marker gene is `GCF_002287175.1_NZ_LMVM01000012.1_cds_WP_069582217.1_1177`, which has genome accession, sequence accession, and the protein name information in it's identifier. 
Once, we have protein and corresponding nucleotide gene sequences from a genome, we used [fetchMG](http://vm-lux.embl.de/~mende/fetchMG/about.html) tool to extract 40 marker gene sequences. Note that fetchMG uses both protein and nucleotide sequences, so keep the name of the two files same, just change the extension to .faa for proteins ana .fna for nucleotide sequences. We used fetchMG v1.0 for our analysis. \
`fetchMG.pl -m extraction -v gene_aa_filename -o output_folder` 

The script get_sequences.py combines the step of downloading genomes and extracting 40 marker genes, just run\
`python get_sequences.py`. \
This should create an output folder with .fna (nucleotide gene sequences) and .faa (protein sequences) files for each of the 40 marker genes. The report.txt has metadata information for the selected gene sequences. Please look at the [output_README.txt](https://github.com/shahnidhi/TIPP_reference_package/blob/master/output_README.txt) file for column headers.  






