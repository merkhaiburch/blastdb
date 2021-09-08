# -------------------------------------------------------------
# Merritt Burch
#
# mbb262@cornell.edu
#
# 2018-11-21
#
# Script to parse maize gff3 file to obtain the most
#   most conserved CDS sequence's 5'UTR-3'UTR start and
#   start-stop of each exon, where overlaps are merged
#
# Input: .gff3 file for maize
#         .txt file containing most conserved maize CDSs
#
# Output: .bed file containing 5' and last 3' UTR positions
#         .bed file containig start and stop pos. of each exon
# -------------------------------------------------------------

# Set working directory
setwd("~/Box Sync/Cornell_PhD/PHG")

# Load in required packages
library(ape)
library(GenomicRanges)

# Pulled data from
# ftp://ftp.ensemblgenomes.org/pub/plants/release-40/gff3/zea_mays

# Load in gene annotation file
gff <- read.gff("Zea_mays.AGPv4.40.chr.gff3")

# Load in list of most conserved CDSs
conserved <- read.table("B73_BLAST_ConservedTranscriptIDs.txt", header = FALSE)


# -----------------------------
# Get ref-ranges for each exon
# -----------------------------

# Get only CDS (each start and stop for each exon)
gff2 <- gff[which(gff$type=="CDS"),]

# Extract only the CDS with transcript IDs (and remove other stuff)
gff2$attributes <- gsub(";Name=.*$", "", gff2$attributes)
gff2$attributes <- gsub(";protein_id=.*$", "", gff2$attributes)
gff2$attributes <- gsub(".*Parent=transcript:", "", gff2$attributes)
gff2$attributes <- gsub(";Parent.*$", "", gff2$attributes)
gff2$attributes <- gsub("ID=CDS:", "", gff2$attributes)

# Subset gff2 exons to only include most conserved CDS
conservedExons = subset(gff2, attributes %in% conserved$V1)

# ---------------------------------------
# Merge overlappping ref-ranges (if-any)
# ---------------------------------------

# Create testing data.frame
exons = conservedExons

# Make genomic range format from data.frame
exons = makeGRangesFromDataFrame(exons,
                                 start.field = "start",
                                 end.field = "end",
                                 seqnames.field = "seqid",
                                 keep.extra.columns = TRUE)

# Reduce ref-ranges (merge overlapping ref-ranges)
exons = reduce(exons)

# Turn this into a data frame and extract relavent columns
exons = as.data.frame(exons)

# Format into .bed format
bed <- exons[,1:3]
bed$start <- as.numeric(bed$start)
bed$start <- bed$start - 1

# Sort data again
bed <- bed[order(bed$seqnames),]

# Write tab delim file
write.table(bed, file = "mbb262_startStop_Exons_refRanges_mergedOverlapping.bed", 
            row.names = FALSE, col.names = FALSE, quote = FALSE,
            sep = "\t")


# -----------------------------
# Get ref-ranges for 5'-3'UTRs
# -----------------------------

# Make a copy of the gff file
UTRs = gff

# Extract only transcript IDs (and remove other stuff) in attributes column
# UTRs$attributes <- gsub(";Name=.*$", "", UTRs$attributes)
# UTRs$attributes <- gsub(";protein_id=.*$", "", UTRs$attributes)
UTRs$attributes <- gsub(".*Parent=transcript:", "", UTRs$attributes) #Parent=transcript:
# UTRs$attributes <- gsub(";Parent.*$", "", UTRs$attributes)
# UTRs$attributes <- gsub("ID=CDS:", "", UTRs$attributes)

# Subset out the 5' and 3' UTR positions
five_utrs <- UTRs[which(UTRs$type=="five_prime_UTR"),]
three_utrs <- UTRs[which(UTRs$type=="three_prime_UTR"),]

# Subset UTRs to only include most conserved CDS
five_utrs <- subset(five_utrs, attributes %in% conserved$V1)
three_utrs <- subset(three_utrs, attributes %in% conserved$V1)

# Take starting position of 5'UTR
five_prime_positiveStrand = five_utrs[which(five_utrs$strand=="+"),]
five_prime_positiveStrand = 
five_prime_negativeStrand = five_utrs[which(five_utrs$strand=="-"),]


# ---------------------------------------
# Merge overlappping ref-ranges (if-any)
# ---------------------------------------

# Create working data.frame
five_utrs_noOverlap = five_utrs
three_utrs_noOverlap = three_utrs

# Make genomic range format from data.frame
five_utrs_noOverlap = makeGRangesFromDataFrame(five_utrs_noOverlap,start.field = "start",end.field = "end",
                                 seqnames.field = "seqid",keep.extra.columns = TRUE)

three_utrs_noOverlap = makeGRangesFromDataFrame(three_utrs_noOverlap, start.field = "start", end.field = "end",
                                 seqnames.field = "seqid",keep.extra.columns = TRUE)

# Reduce ref-ranges (merge overlapping ref-ranges)
five_utrs_noOverlap = reduce(five_utrs_noOverlap)
three_utrs_noOverlap = reduce(three_utrs_noOverlap)

# Turn this into a data frame and extract relavent columns
five_utrs_noOverlap = as.data.frame(five_utrs_noOverlap)
three_utrs_noOverlap = as.data.frame(three_utrs_noOverlap)

# Format into .bed format
bed_five <- five_utrs_noOverlap[,1:3] #5
bed_five$start <- as.numeric(bed_five$start) #5
bed_five$start <- bed_five$start - 1 #5

bed_three <- three_utrs_noOverlap[,1:3] #3
bed_three$start <- as.numeric(bed_three$start) #3
bed_three$start <- bed_three$start - 1 #3

# Sort data again
bed_five <- bed_five[order(bed_five$seqnames),]
bed_three <- bed_three[order(bed_three$seqnames),]

# Write tab delim file
write.table(bed_five, file = "mbb262_fivePrime_refRanges_mergedOverlapping.bed", 
            row.names = FALSE, col.names = FALSE, quote = FALSE,
            sep = "\t")

write.table(bed_three, file = "mbb262_startStop_Exons_refRanges_mergedOverlapping.bed", 
            row.names = FALSE, col.names = FALSE, quote = FALSE,
            sep = "\t")








