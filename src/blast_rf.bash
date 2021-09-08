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
makeblastdb -in Brachypodium_distachyon.Brachypodium_distachyon_v3.0.cds.all.fa -input_type fasta -dbtype nucl -title braRef -out Brachypodium.db
makeblastdb -in Oryza_sativa.IRGSP-1.0.cds.all.fa -input_type fasta -dbtype nucl -title riceRef -out Rice.db
makeblastdb -in Sorghum_bicolor.Sorghum_bicolor_NCBIv3.cds.all.fa -input_type fasta -dbtype nucl -title sorRef -out Sorghum.db

# Combine the databases
blastdb_aliastool -dblist "Brachypodium.db Rice.db Sorghum.db" -dbtype nucl -out three_all -title "Three Genomic" 

# Download maize CDS seqeunces
#wget ftp://ftp.ensemblgenomes.org/pub/plants/release-41/fasta/zea_mays/cds/Zea_mays.B73_RefGen_v4.cds.all.fa.gz
#gzip -d *gz
#mv Zea_mays.B73_RefGen_v4.cds.all.fa maize_cds.fa

# Run BLAST
# Tabular output format
blastn -db ./three_all -query ./Zea_mays.B73_RefGen_v4.cds.all_changedIDs.fa -out matches_tabular_option6_11-7.txt -outfmt "6 qseqid qgi qacc qdescr ssequid sacc sallseqid evalue bitscore length pident nident"
