# ----------------------------------------
# Merritt Burch
#
# 2018-10-31
#
# mbb262@cornell.edu
#
# Script to create custom BLAST database
# containing three genomes
# ---------------------------------------

# Create the databases individually
makeblastdb -in /home/mbb262/blastRefranges/genomes/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.cds.all.fa -input_type fasta -dbtype nucl -title braRef -out Brachypodium.db

makeblastdb -in /home/mbb262/blastRefranges/genomes/Oryza_sativa.IRGSP-1.0.cds.all.fa -input_type fasta -dbtype nucl -title riceRef -out Rice.db

makeblastdb -in /home/mbb262/blastRefranges/genomes/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.cds.all.fa -input_type fasta -dbtype nucl -title sorghumRef -out Sorghum.db

makeblastdb -in /home/mbb262/blastRefranges/genomes/Setaria_italica.Setaria_italica_v2.0.cds.all.fa -input_type fasta -dbtype nucl -title setariaRef -out Setaria.db

# Combine the databases
blastdb_aliastool -dblist "Brachypodium.db Rice.db Sorghum.db Setaria.db" -dbtype nucl -out blast4Genomes -title "blast4Genomes" 

# Download maize CDS seqeunces
# Format CDS files (whitespace sensitieve in further code)
#  ' cds contig:B73_RefGen_v4:' replace with '-'
#  ' cds chromosome:B73_RefGen_v4:' replace with '-'
#wget ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/zea_mays/cds/
#gzip -d *gz
#mv Zea_mays.B73_RefGen_v4.cds.all.fa maize_cds.fa

# Run BLAST
# Tabular output format
blastn -db ./blast4Genomes -query ./Zea_mays.B73_RefGen_v4.42.cds.all.fa -out matches_tabular_option6_2019_02_01.txt -outfmt "6 qseqid qgi qacc qdescr qlen ssequid sacc sallseqid evalue bitscore score length pident nident mismatch gaps"

# Copy all results back to home directory
cp -a /workdir/mbb262/. /home/mbb262/blastRefranges/blast4Genomes/blast4Genomes_CDSv4.42/
