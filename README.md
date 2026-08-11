# 16S Amplicon Analysis: DADA2 and mothur

## Overview

This repository contains my hands-on work exploring two commonly used workflows for 16S rRNA amplicon sequencing analysis: **DADA2** and **mothur**.

This work was conducted as part of my **Research Internship at Nile University**, with a focus on developing practical experience in bioinformatics workflows and building a foundation for **benchmarking bioinformatics tools and pipelines**.

## Objectives

* Explore the main steps of 16S rRNA amplicon analysis
* Gain hands-on experience with DADA2 and mothur
* Understand the differences between ASV- and OTU-based approaches
* Compare the workflows and their analysis steps
* Develop a foundation for evaluating and benchmarking bioinformatics tools

## Workflows

### DADA2

The DADA2 workflow was used to explore an ASV-based approach to 16S rRNA analysis.

Main steps include:

1. Quality assessment
2. Quality filtering and trimming
3. Error model learning
4. Sequence inference
5. Paired-end read merging
6. Sequence table construction
7. Chimera removal
8. Taxonomic assignment
9. Downstream analysis

### mothur

The mothur MiSeq SOP was used to explore an OTU-based approach to 16S rRNA analysis.

Main steps include:

1. Quality assessment
2. Read processing
3. Paired-end read assembly
4. Quality filtering
5. Chimera detection
6. Sequence clustering
7. OTU generation
8. Taxonomic classification
9. Downstream analysis

## DADA2 vs. mothur

One of the main concepts explored in this work is the difference between **ASVs (Amplicon Sequence Variants)** and **OTUs (Operational Taxonomic Units)**.

DADA2 infers exact sequence variants using an error model, while traditional OTU-based approaches cluster sequences according to a similarity threshold.

Therefore, the resulting features have different biological definitions and should not be treated as directly equivalent.

## Key Learning Outcomes

Through this work, I gained practical experience with:

* 16S rRNA amplicon sequencing
* Quality control and preprocessing
* ASV-based analysis
* OTU-based analysis
* Taxonomic assignment
* Microbial community analysis
* Comparing bioinformatics workflows

## Research Context

This work is part of my ongoing research internship at **Nile University**, where I am developing my understanding of **bioinformatics benchmarking** and learning how different tools and workflows can be systematically evaluated and compared.

## References

* DADA2 Tutorial: https://benjjneb.github.io/dada2/tutorial.html
* mothur MiSeq SOP: https://mothur.org/wiki/miseq_sop/
