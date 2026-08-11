# mothur 16S rRNA Amplicon Analysis

## Overview

This folder contains my hands-on implementation of the **mothur MiSeq SOP** for 16S rRNA amplicon sequencing analysis.

This work was conducted as part of my **Research Internship at Nile University**, where I am developing practical experience in bioinformatics workflows and building a foundation for **bioinformatics benchmarking**.

The mothur workflow processes paired-end 16S rRNA amplicon reads through quality control, read assembly, alignment, chimera detection, taxonomic classification, and OTU generation.

## Objective

The main objectives of this analysis were to:

* Explore the mothur MiSeq SOP for 16S rRNA amplicon data
* Gain hands-on experience with mothur
* Understand paired-end read processing and quality filtering
* Perform sequence alignment and filtering
* Detect and remove chimeric sequences
* Assign taxonomy to sequences
* Generate OTUs and an OTU abundance table
* Evaluate the workflow using a Mock community
* Understand the differences between OTU-based and ASV-based approaches

## Workflow

The analysis follows the main steps of the mothur MiSeq SOP:

```text
Raw Paired-End FASTQ Reads
          ↓
Create File
          ↓
Make Contigs
          ↓
Quality Filtering
          ↓
Dereplication
          ↓
Reference Alignment
          ↓
Alignment Filtering
          ↓
Pre-clustering
          ↓
Chimera Detection
          ↓
Taxonomic Classification
          ↓
Remove Unwanted Lineages
          ↓
Mock Community Evaluation
          ↓
Remove Mock Sample
          ↓
OTU Clustering
          ↓
OTU Abundance Table
          ↓
Downstream Analysis
```

## Main Analysis Steps

### 1. Create the Input File

The `make.file` command is used to create a file describing the paired-end FASTQ files and their corresponding samples.

```bash
mothur "#make.file(inputdir=., type=fastq, prefix=stability)"
```

This provides mothur with the information needed to identify forward and reverse reads for each sample.

### 2. Assemble Paired-End Reads

Forward and reverse reads are assembled into contigs using:

```bash
mothur "#make.contigs(file=stability.files)"
```

This step combines the paired-end reads based on their overlapping region.

### 3. Quality Filtering

The generated contigs are inspected and filtered based on sequence quality.

Sequences containing ambiguous bases, excessive length, or long homopolymers are removed.

The main filtering criteria used in the workflow include:

* `maxambig=0`
* `maxlength=275`
* `maxhomop=8`

### 4. Dereplication

Identical sequences are collapsed using:

```bash
mothur "#unique.seqs()"
```

This reduces redundant sequence processing while preserving abundance information through the count table.

### 5. Reference Alignment

Sequences are aligned against a customized **SILVA reference alignment**.

The reference is first restricted to the target region and then used to align the processed sequences.

### 6. Alignment Filtering

Sequences are screened to retain reads covering the expected region of the alignment.

Alignment columns containing only gaps are also removed.

This helps ensure that sequences being compared cover a consistent region.

### 7. Pre-clustering

The `pre.cluster` command is used to merge sequences that differ by a small number of nucleotides.

In this workflow:

```bash
diffs=2
```

is used to reduce the effect of sequencing and PCR errors before chimera detection.

### 8. Chimera Detection

Potential chimeric sequences are identified using:

```bash
chimera.vsearch
```

Chimeric sequences are removed before taxonomic classification and OTU generation.

### 9. Taxonomic Classification

Sequences are classified using a mothur-formatted reference database and taxonomy file.

The workflow uses:

```bash
classify.seqs
```

to assign taxonomy to the processed sequences.

### 10. Remove Unwanted Lineages

Sequences classified as unwanted groups are removed before downstream analysis.

The workflow removes:

* Chloroplast
* Mitochondria
* Unknown
* Archaea
* Eukaryota

This step focuses the analysis on the bacterial sequences of interest.

### 11. Mock Community Evaluation

A Mock community is used to evaluate the sequencing and processing workflow.

The `seq.error` command compares observed sequences with expected reference sequences to estimate sequencing error.

This provides a way to assess how well the workflow recovers the expected Mock community sequences.

### 12. Remove the Mock Sample

After evaluation, the Mock sample is removed before downstream microbial community analysis.

This prevents the artificial Mock community from being included in analyses of the biological samples.

### 13. OTU Clustering

Sequences are clustered into **Operational Taxonomic Units (OTUs)**.

A 3% distance cutoff is used:

```bash
mothur "#dist.seqs(fasta=final.fasta, cutoff=0.03)"
```

followed by clustering.

A 0.03 distance corresponds to the traditional **97% sequence similarity** OTU definition.

### 14. Generate the OTU Table

The `make.shared` command is used to generate an OTU abundance table containing the number of sequences assigned to each OTU in each sample.

The resulting shared file can be used for downstream microbial community analysis.

## OTU-Based Approach

A key concept explored in this workflow is the use of **Operational Taxonomic Units (OTUs)**.

In a traditional OTU-based workflow, sequences are grouped according to a sequence similarity threshold.

In this analysis, sequences are clustered at a 3% distance threshold, corresponding to approximately 97% similarity.

This differs from the **ASV-based approach used by DADA2**, which attempts to infer exact sequence variants by modeling sequencing errors.

Therefore, OTUs and ASVs represent **different biological feature definitions** and should not be treated as directly equivalent.

## Key Learning Outcomes

Through this workflow, I gained hands-on experience with:

* 16S rRNA amplicon sequencing
* Paired-end read processing
* Quality filtering
* Sequence dereplication
* Reference sequence alignment
* Alignment filtering
* Pre-clustering
* Chimera detection
* Taxonomic classification
* Mock community evaluation
* OTU clustering
* OTU abundance table generation
* Comparing OTU-based and ASV-based approaches

## Tools and References

| Tool / Resource | Purpose                                                       |
| --------------- | ------------------------------------------------------------- |
| mothur          | 16S rRNA sequence processing and microbial community analysis |
| SILVA           | Reference sequence alignment and taxonomy                     |
| VSEARCH         | Chimera detection                                             |
| RDP             | Taxonomic classification reference                            |
| Mock community  | Workflow evaluation                                           |

## Repository Files

```text
mothur/
│
├── README.md
│
├── scripts/
│   └── mothur_analysis.sh
```

## Research Context

This work is part of my ongoing **Research Internship at Nile University**, where I am developing my understanding of **bioinformatics benchmarking** and learning how different tools and workflows can be systematically evaluated and compared.

The mothur workflow provides an important comparison point with the DADA2 workflow in this repository, particularly because of the difference between **OTU-based clustering and ASV-based inference**.

## Reference

The workflow is based on the official mothur MiSeq SOP:

* [mothur MiSeq SOP](https://mothur.org/wiki/miseq_sop/)
