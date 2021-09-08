# ----------------------
# Merritt Burch
# 2019-02-11
# mbb262@cornell.edu
# ----------------------

# ---------------------------------
# BLASTP of maize proteins against
# Sorghum, seteria, brachy, rice
# ---------------------------------

## -----------------
## Get protein data
## -----------------

#wget ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/zea_mays/pep/Zea_mays.B73_RefGen_v4.pep.all.fa.gz

wget ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/brachypodium_distachyon/pep/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.pep.all.fa.gz

wget ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/oryza_sativa/pep/Oryza_sativa.IRGSP-1.0.pep.all.fa.gz

wget ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/sorghum_bicolor/pep/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.pep.all.fa.gz

wget ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/setaria_italica/pep/Setaria_italica.Setaria_italica_v2.0.pep.all.fa.gz


## -----------------
## Unzip all files
## -----------------
gzip -d *.gz


## -----------------
## Make databases
## -----------------
makeblastdb -in ./Brachypodium_distachyon.Brachypodium_distachyon_v3.0.pep.all.fa -input_type fasta -dbtype prot -title braRefprot -out BrachypodiumProt.db

makeblastdb -in ./Oryza_sativa.IRGSP-1.0.pep.all.fa -input_type fasta -dbtype prot -title riceRefprot -out RiceProt.db

makeblastdb -in ./Sorghum_bicolor.Sorghum_bicolor_NCBIv3.pep.all.fa -input_type fasta -dbtype prot -title sorghumRefprot -out SorghumProt.db

makeblastdb -in ./Setaria_italica.Setaria_italica_v2.0.pep.all.fa -input_type fasta -dbtype prot -title setariaRefprot -out SetariaProt.db



# Format maize CDS files (whitespace sensitieve in further code)
#  ' pep contig:B73_RefGen_v4:' replace with '-'
#  ' pep chromosome:B73_RefGen_v4:' replace with '-'
#wget ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/zea_mays/cds/
#gzip -d *gz
#mv Zea_mays.B73_RefGen_v4.cds.all.fa maize_cds.fa


# Combine the databases
blastdb_aliastool -dblist "BrachypodiumProt.db RiceProt.db SorghumProt.db SetariaProt.db" -dbtype prot -out blast4GenomesProt -title "blast4GenomesProtein"

## -----------------
## Run BLAST
## -----------------

# Tabular output format
blastp -db ./blast4GenomesProt -query ./Zea_mays.B73_RefGen_v4.42.pep.all.fa -out matches_tabular_option6_2019_02_11_protein.txt -outfmt "6 qseqid qgi qacc qdescr qlen ssequid sacc sallseqid evalue bitscore score length pident nident mismatch gaps"

# Copy all results back to home directory
cp -a /workdir/mbb262/. /home/mbb262/blastRefranges/blast4Genomes/blast4Genomes_CDSv4.42_protein/

