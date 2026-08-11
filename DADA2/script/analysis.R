# ============================================================
# DADA2 16S rRNA Amplicon Analysis
# ============================================================

# 1. Load DADA2

library(dada2)


# ============================================================
# 2. Define the Input Data Directory
# ============================================================

# Set the path containing the raw paired-end FASTQ files
path <- "data"

# List all files in the input directory
list.files(path)


# ============================================================
# 3. Identify Forward and Reverse Reads
# ============================================================

# Identify forward reads (R1) and reverse reads (R2)
# and sort them so that corresponding pairs are in the same order
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))


# ============================================================
# 4. Extract Sample Names
# ============================================================

# Extract sample names from the FASTQ filenames
# The first part of each filename is used as the sample name
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Display the sample names
sample.names


# ============================================================
# 5. Assess Read Quality
# ============================================================

# Visualize quality profiles of the first two forward reads
plotQualityProfile(fnFs[1:2])

# Visualize quality profiles of the first two reverse reads
plotQualityProfile(fnRs[1:2])


# ============================================================
# 6. Define Output Paths for Filtered Reads
# ============================================================

# Create file paths for the filtered forward and reverse reads
filtFs <- file.path(path, "filtered",
                    paste0(sample.names, "_F_filt.fastq.gz"))

filtRs <- file.path(path, "filtered",
                    paste0(sample.names, "_R_filt.fastq.gz"))

# Assign sample names to the filtered read file paths
names(filtFs) <- sample.names
names(filtRs) <- sample.names


# ============================================================
# 7. Filter and Trim Reads
# ============================================================

# Filter and trim the raw reads based on quality criteria
#
# truncLen:
#   Truncate forward reads at 240 bp and reverse reads at 160 bp
#
# maxN:
#   Allow no ambiguous bases (N)
#
# maxEE:
#   Allow a maximum expected error of 2 for each read
#
# truncQ:
#   Truncate reads when a quality score of 2 or less is encountered
#
# rm.phix:
#   Remove reads matching PhiX sequences
#
# compress:
#   Compress the output FASTQ files
#
# multithread:
#   Use multiple CPU threads for faster processing
#   On Windows, set multithread=FALSE if required

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(240, 160),
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

# Display the first rows of the filtering summary
head(out)


# ============================================================
# 8. Learn the Error Rates
# ============================================================

# Learn the sequencing error model from the filtered
# forward and reverse reads
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Visualize the learned error models
plotErrors(errF, nominalQ = TRUE)


# ============================================================
# 9. Infer Amplicon Sequence Variants (ASVs)
# ============================================================

# Apply the learned error models to infer sequence variants
# separately for forward and reverse reads
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)


# ============================================================
# 10. Merge Paired-End Reads
# ============================================================

# Merge the denoised forward and reverse reads
# based on their overlapping regions
mergers <- mergePairs(
  dadaFs, filtFs,
  dadaRs, filtRs,
  verbose = TRUE
)

# Inspect the merger results for the first sample
head(mergers[[1]])


# ============================================================
# 11. Construct the Sequence Table
# ============================================================

# Create a sequence table containing the inferred
# ASVs and their abundances across samples
seqtab <- makeSequenceTable(mergers)

# Display the dimensions of the sequence table
dim(seqtab)


# ============================================================
# 12. Inspect Sequence Length Distribution
# ============================================================

# Examine the distribution of the lengths of inferred sequences
table(nchar(getSequences(seqtab)))


# ============================================================
# 13. Remove Chimeric Sequences
# ============================================================

# Identify and remove chimeric sequences using the consensus method
seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

# Display the dimensions after chimera removal
dim(seqtab.nochim)

# Calculate the proportion of reads remaining after
# chimera removal
sum(seqtab.nochim) / sum(seqtab)


# ============================================================
# 14. Track Reads Through the Pipeline
# ============================================================

# Function to calculate the number of unique reads
getN <- function(x) sum(getUniques(x))

# Combine read counts from each major processing step
# to track how many reads remain throughout the workflow
track <- cbind(
  out,
  sapply(dadaFs, getN),
  sapply(dadaRs, getN),
  sapply(mergers, getN),
  rowSums(seqtab.nochim)
)

# If processing a single sample, remove the sapply calls:
# e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

# Assign descriptive column names to the tracking table
colnames(track) <- c(
  "input",
  "filtered",
  "denoisedF",
  "denoisedR",
  "merged",
  "nonchim"
)

# Assign sample names as row names
rownames(track) <- sample.names

# Display the first rows of the tracking table
head(track)


# ============================================================
# 15. Assign Taxonomy
# ============================================================

# Assign taxonomy to the non-chimeric ASVs using
# the SILVA reference database
taxa <- assignTaxonomy(
  seqtab.nochim,
  "silva_nr99_v138.2_toGenus_trainset.fa.gz",
  multithread = TRUE
)

# Create a copy of the taxonomy table for display
# without sequence identifiers as row names
taxa.print <- taxa
rownames(taxa.print) <- NULL

# Display the first rows of the taxonomy table
head(taxa.print)


# ============================================================
# 16. Evaluate the Mock Community
# ============================================================

# Extract ASVs detected in the Mock community sample
unqs.mock <- seqtab.nochim["Mock",]

# Keep only ASVs that are present in the Mock sample
# and sort them by abundance
unqs.mock <- sort(
  unqs.mock[unqs.mock > 0],
  decreasing = TRUE
)

# Report the number of ASVs detected in the Mock community
cat(
  "DADA2 inferred",
  length(unqs.mock),
  "sample sequences present in the Mock community.\n"
)

# Load the expected reference sequences for the Mock community
mock.ref <- getSequences(
  file.path(path, "HMP_MOCK.v35.fasta")
)

# Check how many inferred sequences exactly match
# the expected reference sequences
match.ref <- sum(
  sapply(
    names(unqs.mock),
    function(x) any(grepl(x, mock.ref))
  )
)

# Report the number of exact matches
cat(
  "Of those,",
  sum(match.ref),
  "were exact matches to the expected reference sequences.\n"
)


# ============================================================
# 17. Load Packages for Downstream Analysis
# ============================================================

library(phyloseq)
library(Biostrings)
library(ggplot2)

# Set the default ggplot2 theme
theme_set(theme_bw())


# ============================================================
# 18. Prepare Sample Metadata
# ============================================================

# Extract sample names from the non-chimeric sequence table
samples.out <- rownames(seqtab.nochim)

# Extract subject identifiers from the sample names
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)

# Extract gender information from the subject identifier
gender <- substr(subject, 1, 1)

# Remove the gender character from the subject identifier
subject <- substr(subject, 2, 999)

# Extract the sampling day from the sample name
day <- as.integer(
  sapply(strsplit(samples.out, "D"), `[`, 2)
)

# Create a sample metadata data frame
samdf <- data.frame(
  Subject = subject,
  Gender = gender,
  Day = day
)

# Categorize samples into Early and Late time points
samdf$When <- "Early"
samdf$When[samdf$Day > 100] <- "Late"

# Set sample names as row names
rownames(samdf) <- samples.out


# ============================================================
# 19. Create a Phyloseq Object
# ============================================================

# Combine the ASV table, sample metadata, and taxonomy
# into a phyloseq object for downstream analysis
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  sample_data(samdf),
  tax_table(taxa)
)

# Remove the Mock sample before downstream community analysis
ps <- prune_samples(
  sample_names(ps) != "Mock",
  ps
)


# ============================================================
# 20. Add ASV Sequences to the Phyloseq Object
# ============================================================

# Convert ASV sequence names into DNAStringSet objects
dna <- Biostrings::DNAStringSet(taxa_names(ps))

# Assign the ASV sequences as names
names(dna) <- taxa_names(ps)

# Add the DNA sequences to the phyloseq object
ps <- merge_phyloseq(ps, dna)

# Rename ASVs using simple identifiers
# (ASV1, ASV2, ASV3, ...)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# Display the phyloseq object
ps


# ============================================================
# 21. Alpha Diversity Analysis
# ============================================================

# Calculate and visualize Shannon and Simpson diversity
# across sampling days
plot_richness(
  ps,
  x = "Day",
  measures = c("Shannon", "Simpson"),
  color = "When"
)


# ============================================================
# 22. Prepare Data for Bray-Curtis Analysis
# ============================================================

# Transform abundance data into relative proportions
# before calculating Bray-Curtis distances
ps.prop <- transform_sample_counts(
  ps,
  function(otu) otu / sum(otu)
)


# ============================================================
# 23. NMDS Ordination
# ============================================================

# Perform NMDS ordination using Bray-Curtis distance
ord.nmds.bray <- ordinate(
  ps.prop,
  method = "NMDS",
  distance = "bray"
)

# Visualize the NMDS ordination
plot_ordination(
  ps.prop,
  ord.nmds.bray,
  color = "When",
  title = "Bray NMDS"
)


# ============================================================
# 24. Visualize the Most Abundant Taxa
# ============================================================

# Identify the 20 most abundant taxa
top20 <- names(
  sort(
    taxa_sums(ps),
    decreasing = TRUE
  )
)[1:20]

# Convert abundance values to relative proportions
ps.top20 <- transform_sample_counts(
  ps,
  function(OTU) OTU / sum(OTU)
)

# Keep only the top 20 most abundant taxa
ps.top20 <- prune_taxa(top20, ps.top20)

# Create a stacked bar plot showing taxonomic composition
# across sampling days and time points
plot_bar(
  ps.top20,
  x = "Day",
  fill = "Family"
) +
  facet_wrap(
    ~When,
    scales = "free_x"
  )

