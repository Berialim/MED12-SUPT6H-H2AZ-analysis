# Inputs and execution

The archive is organized by source path to preserve relative imports and directory relationships. Most scripts require manual path configuration. Work on copies of the input data and inspect each script before execution.

| Module | Required inputs |
| --- | --- |
| Proteomics | Original protein differential tables, including `volcannodata.csv` for the RPB1 script; sample intensities/design for primary analysis |
| ChIP-seq | FASTQ/BAM, input controls, hg38 and dm6 indices, exact spike-in counts; peaks and compatible annotations for ChIPseeker |
| PRO-seq | Strand-specific normalized BigWigs, BED12 gene coordinates, gene IDs; prepared `*_pause` TSV columns for plotting |
| TT-seq | `featureCounts.txt`, `database.csv`, `scalingfactors.txt`, orthology table and gene-level results; BED column 4 gene ID and column 7 H2A.Z density for the tertile plot |
| MNase | `med12tss_to_plus1_peak.bed`, `spt6tss_to_plus1_peak.bed`, matching BED6 gene IDs; upstream peak caller still required |
| Mouse RNA-seq | Count matrices and sample metadata for both overexpression cohorts, mm39 reference and annotation |
| Composite score | `h2az/mousetss500_h2az.txt`, `h2az/refGene_mm39.txt.gz`, `med12_oe/results.txt`, `spt6_oe/results.txt`, versioned species-appropriate GMT files |
| External comparison | Published MED12 loss-of-function results and the exact metadata/contrast; verify GSE229062 correspondence before use |

Python dependencies are listed in `requirements-analysis.txt`; these are observed dependencies, not a tested lockfile. R dependencies observed include DESeq2, tidyverse, ggpubr, ggrepel, pheatmap, VennDiagram, clusterProfiler, enrichplot, RColorBrewer and ChIPseeker with hg38 TxDb/organism annotation packages. Other modules may require further dependencies. Shell scripts need a Linux/WSL environment with their referenced tools.

Versions reported in the manuscript (not verified installed): Trim Galore 0.6.10, Bowtie2 2.5.1, Picard 2.27.5, MACS2 2.2.7/2.2.7.1, deepTools 3.5.1, STAR 2.7.10a, SAMtools 1.17, featureCounts 2.0.3, DESeq2 1.40.0 and GSEApy 1.0.4. Record actual versions and references for every future run. The manuscript's human effective genome size must not be applied automatically to mouse.

Before inference, verify biological replicates, sample identity, reference assemblies, strand orientation, gene-ID uniqueness, normalization, contrast direction and multiple-testing procedures. No sample design or numerical result is invented in this repository.
