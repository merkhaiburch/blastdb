# ------------------------------------------------------
# Merritt Burch
#
# mbb262@cornell.edu
#
# 2018-11-07
#
# Script to parse three genome BLAST results
# for summary statstics, 
# This version has updated gene IDs that include
# the start and stop positions
#
# Input: blast tabular output from alignment of 
#        B73 CDSs with Brachy, sorghum, and rice CDSs
#
# Output: .txt file containing most conserved CDSs
#         .bed file containig testable ref-ranges for PHG
#
# Flavor: Positions of ATG-start codon and terminator with 
#         merged overlaps
# --------------------------------------------------------

# Set working directory
setwd("~/Box Sync/Cornell_PhD/PHG")

# Load file containg positions of start codon - stop codon
blast <- read.table("matches_tabular_option6_11-8.txt", header=FALSE)


# ------------------
#  Parse out file
# ------------------

#Separate first column into general gene, chromosome, start,stop
# Not pretty but this will work
blast$gene = blast$V1
blast$chr = blast$V1
blast$start = blast$V1
blast$stop = blast$V1

# Format
blast$gene = gsub("-[0-9]:.*", "", blast$gene)
blast$gene = gsub("-[0-9][0-9]:.*", "", blast$gene)
#blast$gene = gsub("-B73V4_.*", "", blast$gene)

# Remove rows of mitochondrial and chloroplast genes
blast <- blast[!grepl("-Mt:", blast$gene),]
blast <- blast[!grepl("-Pt:", blast$gene),]
blast <- blast[!grepl("B73V4_ctg", blast$gene),]

# More general formatting
blast$transcript = blast$gene
blast$gene = gsub("_T[0-9][0-9][0-9]", "", blast$gene)
blast$gene = gsub("_T[0-9][0-9]", "", blast$gene)
blast$chr = gsub(".*[0-9]-", "", blast$chr)
blast$chr = gsub(":.*", "", blast$chr)

# Format start position
#blast$start = gsub(".*_ctg[0-9][0-9][0-9]:", "", blast$start)
#blast$start = gsub(".*_ctg[0-9][0-9]:", "", blast$start)
#blast$start = gsub(".*_ctg[0-9]:", "", blast$start)
blast$start = gsub(".*-[0-9]:", "", blast$start)
blast$start = gsub(".*-10:", "", blast$start)
blast$start = gsub(":.*", "", blast$start)

# Format stop position
blast$stop = gsub(":-1$", "", blast$stop)
blast$stop = gsub(":1$", "", blast$stop)
blast$stop = gsub(".*:", "", blast$stop)

# Rearrange
blast <- blast[,c(11,15,12,13,14,6,7,8,9,10,4)]

# Rename columns
colnames(blast) <- c("gene", "trans", "chr", "five_prime_UTR_start", "three_prime_UTR_stop",
                      "evalue", "bitScore", "rawValue", 
                      "percentMatch", "other", "proteinMatch1")

# -------------------
# Filter out results
# -------------------

# Take top BLAST hit with highest bit score & percent match, lowest evalue
# Subsets by general gene
library(dplyr)
refRangegenes <- blast %>% 
  group_by(gene) %>% 
  filter(evalue==min(evalue) & bitScore==max(bitScore))

# Sort data by gene ID and transcript ID
refRangegenes <- refRangegenes[order(refRangegenes$gene, refRangegenes$trans),]

# Subset out only unique genes (because many transcripts look
#   to be duplicated)
refRangegenes2 <- refRangegenes[!duplicated(refRangegenes$gene),]

# Export only transcript IDs
transcIds <- as.data.frame(refRangegenes2$trans)
write.table(transcIds, file = "B73_BLAST_ConservedTranscriptIDs.txt",
            row.names = FALSE, col.names = FALSE)


# --------------------------------
# Look for overlapping ref-ranges
# --------------------------------

# Create testing data.frame
ranges = refRangegenes2

# Use genomic ranges to merge overlapping ranges
library(GenomicRanges)

# Make genomic range format from data.frame
refRangegenes2 = makeGRangesFromDataFrame(refRangegenes2,
                               start.field = "five_prime_UTR_start",
                               end.field = "three_prime_UTR_stop",
                               seqnames.field = "chr",
                               keep.extra.columns = TRUE)

# Reduce ref-ranges (merge overlapping ref-ranges)
refRangegenes2= reduce(refRangegenes2)

# Turn this into a data frame and extract relavent columns
refRangegenes2 = as.data.frame(refRangegenes2)


# ------------------------------------------------------
# Get CDS start/stop positions for most conserved genes
# ------------------------------------------------------

# Format to Lynn's .bed format
bed <- refRangegenes2[,1:3]
bed$start <- as.numeric(bed$start)
bed$start <- bed$start - 1

# Sort data again
bed <- bed[order(bed$chr),]

# Write tab delim file
write.table(bed, file = "mbb262_UTRsIncluded_refRanges_mergedOverlapping.bed", 
            row.names = FALSE, col.names = FALSE, quote = FALSE,
            sep = "\t")


