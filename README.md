# blastdb
Create a custom blastdb and search sequences against it

Conserved CDSs from B73_RefGen_v4
Content
src/: R code for subsetting BLASTdb hits and creating testable ref-ranges; bash scripts for running BLASTdb

Methods
A BLAST database was created using CDS sequences taken from Oryza sativa (IRGSP-1.0), Sorghum bicolor (NCBIv3), and Brachypodium distachyon (v3.0). B73 CDS sequences were searched against this BLAST DB, top hits were parsed using R. Non-overlapping start and stop positions for each exon, the 5' UTR and 3' UTR, and the transcriptional start site and the transcriptional stop sites are in tab separated bed files for testing in the PHG.

Contacts
Merritt Burch (mbb262@cornell.edu)