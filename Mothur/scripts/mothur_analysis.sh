#!/bin/bash

# ============================================================
# mothur MiSeq SOP - 16S rRNA Amplicon Analysis
# ============================================================
#
# Description:
# This script follows the main workflow described in the
# mothur MiSeq SOP for paired-end 16S rRNA amplicon data.
#
# NOTE:
# The reference files, input data, sample names, and some
# parameters may need to be adapted to your dataset.
# ============================================================


# ============================================================
# 1. Create the mothur input file
# ============================================================

mothur "#make.file(inputdir=., type=fastq, prefix=stability)"


# ============================================================
# 2. Assemble Paired-End Reads
# ============================================================

mothur "#make.contigs(file=stability.files)"


# ============================================================
# 3. Inspect Sequence Quality
# ============================================================

mothur "#summary.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table)"


# ============================================================
# 4. Quality Filtering
# ============================================================

mothur "#screen.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table, maxambig=0, maxlength=275, maxhomop=8)"


# ============================================================
# 5. Remove Duplicate Sequences
# ============================================================

mothur "#unique.seqs(fasta=stability.trim.contigs.good.fasta, count=stability.contigs.good.count_table)"


# ============================================================
# 6. Prepare the Reference Alignment
# ============================================================

mothur "#pcr.seqs(fasta=silva.bacteria.fasta, start=11895, end=25318, keepdots=F)"

mothur "#rename.file(input=silva.bacteria.pcr.fasta, new=silva.v4.fasta)"

mothur "#summary.seqs(fasta=silva.v4.fasta)"


# ============================================================
# 7. Align Sequences to the Reference
# ============================================================

mothur "#align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=silva.v4.fasta)"

mothur "#summary.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.contigs.good.count_table)"


# ============================================================
# 8. Keep Sequences Covering the Same Region
# ============================================================

mothur "#screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.contigs.good.count_table, start=1968, end=11550)"

mothur "#summary.seqs(fasta=current, count=current)"


# ============================================================
# 9. Filter the Alignment
# ============================================================

mothur "#filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=.)"


# ============================================================
# 10. Remove Duplicate Sequences Again
# ============================================================

mothur "#unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.contigs.good.unique.good.filter.count_table)"


# ============================================================
# 11. Pre-cluster Sequences
# ============================================================

mothur "#pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table, diffs=2)"


# ============================================================
# 12. Detect and Remove Chimeras
# ============================================================

mothur "#chimera.vsearch(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=T)"

mothur "#summary.seqs(fasta=current, count=current)"


# ============================================================
# 13. Assign Taxonomy
# ============================================================

mothur "#classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.count_table, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax)"


# ============================================================
# 14. Remove Unwanted Lineages
# ============================================================

mothur "#remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)"

mothur "#summary.tax(taxonomy=current, count=current)"


# ============================================================
# 15. Evaluate the Mock Community
# ============================================================

mothur "#get.groups(count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.fasta, groups=Mock)"

mothur "#seq.error(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, reference=HMP_MOCK.v35.fasta, aligned=F)"


# ============================================================
# 16. Remove the Mock Sample
# ============================================================

mothur "#remove.groups(count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.fasta, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pds.wang.pick.taxonomy, groups=Mock)"

mothur "#rename.file(fasta=current, count=current, taxonomy=current, prefix=final)"


# ============================================================
# 17. Generate OTUs
# ============================================================

mothur "#dist.seqs(fasta=final.fasta, cutoff=0.03)"

mothur "#cluster(column=final.dist, count=final.count_table)"


# ============================================================
# 18. Create the OTU Abundance Table
# ============================================================

mothur "#make.shared(list=final.opti_mcc.list, count=final.count_table, label=0.03)"


# ============================================================
# 19. Assign Consensus Taxonomy to OTUs
# ============================================================

mothur "#classify.otu(list=final.opti_mcc.list, count=final.count_table, taxonomy=final.taxonomy, label=0.03)"


# ============================================================
# 20. Count Sequences per Sample
# ============================================================

mothur "#count.groups(shared=final.opti_mcc.shared)"


# ============================================================
# 21. Subsample the Dataset
# ============================================================
#
# NOTE:
# 2403 is the value used in the tutorial dataset.
# Replace it with the appropriate minimum sample size
# for your own dataset.

mothur "#sub.sample(shared=final.opti_mcc.shared, size=2403)"


# ============================================================
# 22. Rarefaction Analysis
# ============================================================

mothur "#rarefaction.single(shared=final.opti_mcc.shared, calc=sobs, freq=100)"


# ============================================================
# End of Workflow
# ============================================================

echo "mothur 16S rRNA analysis workflow completed."
