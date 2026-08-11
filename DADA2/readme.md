# DADA2 16S rRNA Amplicon Analysis

## Overview

This folder contains my hands-on implementation of a **DADA2-based workflow for 16S rRNA amplicon sequencing analysis**.

This work was conducted as part of my **Research Internship at Nile University**, where I am developing practical experience in bioinformatics workflows and building a foundation for **bioinformatics benchmarking**.

DADA2 is an ASV-based approach that models sequencing errors to infer exact **Amplicon Sequence Variants (ASVs)** from amplicon sequencing data.

## Objective

The main objectives of this analysis were to:

* Explore the DADA2 workflow for 16S rRNA amplicon data
* Perform quality assessment and read preprocessing
* Learn how DADA2 infers ASVs
* Perform paired-end read merging
* Remove chimeric sequences
* Assign taxonomy using a reference database
* Evaluate the results using a Mock community
* Perform downstream microbial community analysis using phyloseq

## Workflow

The analysis follows these main steps:

```text
Raw FASTQ Reads
      ↓
Quality Assessment
      ↓
Filtering and Trimming
      ↓
Learn Error Rates
      ↓
Denoising / ASV Inference
      ↓
Paired-End Read Merging
      ↓
Sequence Table
      ↓
Chimera Removal
      ↓
Taxonomic Assignment
      ↓
Mock Community Evaluation
      ↓
Phyloseq Analysis
      ↓
Diversity & Community Analysis
```

## Main Analysis Steps

### 1. Quality Assessment

The quality profiles of forward and reverse reads were inspected using:

```r
plotQualityProfile()
```

This helps determine appropriate trimming and filtering parameters.

### 2. Filtering and Trimming

Reads were filtered using:

```r
filterAndTrim()
```

The filtering criteria included:

* Maximum expected errors: 2
* No ambiguous bases
* Quality-based truncation
* PhiX removal

Forward and reverse reads were truncated according to their quality profiles.

### 3. Error Model Learning

DADA2 learns sequencing error rates from the filtered reads using:

```r
learnErrors()
```

Separate error models were generated for forward and reverse reads.

### 4. ASV Inference

Denoising and ASV inference were performed using:

```r
dada()
```

This step models sequencing errors and distinguishes true biological sequence variants from sequencing errors.

### 5. Paired-End Read Merging

Forward and reverse reads were merged using:

```r
mergePairs()
```

The overlapping regions of the paired-end reads are used to reconstruct the full amplicon sequence.

### 6. Sequence Table Construction

The inferred sequences were organized into an ASV abundance table using:

```r
makeSequenceTable()
```

Rows represent samples, while columns represent inferred ASVs.

### 7. Chimera Removal

Chimeric sequences were identified and removed using:

```r
removeBimeraDenovo()
```

This helps reduce artificial sequences generated during PCR amplification.

### 8. Read Tracking

Reads were tracked throughout the workflow to monitor how many sequences remained after each processing step.

The main stages tracked were:

| Stage     | Description                             |
| --------- | --------------------------------------- |
| Input     | Raw reads                               |
| Filtered  | Reads remaining after quality filtering |
| DenoisedF | Denoised forward reads                  |
| DenoisedR | Denoised reverse reads                  |
| Merged    | Successfully merged paired reads        |
| Nonchim   | Reads remaining after chimera removal   |

### 9. Taxonomic Assignment

Taxonomy was assigned to the non-chimeric ASVs using a **SILVA reference database** with:

```r
assignTaxonomy()
```

### 10. Mock Community Evaluation

A Mock community was used to evaluate the inferred ASVs.

The inferred sequences were compared against the expected reference sequences to determine how many DADA2-inferred sequences were exact matches to the expected Mock community sequences.

## Downstream Analysis

The resulting ASV table was imported into **phyloseq** for downstream microbial community analysis.

The analysis included:

### Alpha Diversity

Shannon and Simpson diversity indices were calculated to explore within-sample diversity.

### Beta Diversity

Relative abundance data were used to calculate **Bray-Curtis distances**, followed by **NMDS ordination** to explore differences in microbial community composition.

### Taxonomic Composition

The most abundant taxa were visualized using a stacked bar plot at the **Family** level.

## Tools and Packages

| Tool / Package | Purpose                                        |
| -------------- | ---------------------------------------------- |
| DADA2          | 16S rRNA sequence processing and ASV inference |
| phyloseq       | Microbial community analysis                   |
| Biostrings     | Biological sequence manipulation               |
| ggplot2        | Data visualization                             |
| SILVA          | Taxonomic reference database                   |

## Key Learning Outcomes

Through this workflow, I gained hands-on experience with:

* 16S rRNA amplicon sequencing
* FASTQ quality assessment
* Quality filtering and trimming
* Sequencing error modeling
* ASV inference
* Paired-end read merging
* Chimera detection
* Taxonomic assignment
* Mock community evaluation
* Alpha and beta diversity analysis
* Phyloseq-based downstream analysis

## ASV-Based Approach

A key concept explored in this workflow is the use of **Amplicon Sequence Variants (ASVs)**.

Unlike traditional OTU clustering, DADA2 attempts to resolve exact biological sequence variants by modeling and correcting sequencing errors.

This makes ASVs useful for comparing sequence-level variation across samples and provides an important contrast with the **OTU-based approach used in mothur**, which is explored in the corresponding folder of this repository.

## Repository Files

```text
DADA2/
│
├── README.md
├── scripts/
│   └── dada2_analysis.R
└── figures/
```

## Reference

The workflow was based on the DADA2 tutorial:

* [DADA2 MiSeq Tutorial](https://benjjneb.github.io/dada2/tutorial.html)
