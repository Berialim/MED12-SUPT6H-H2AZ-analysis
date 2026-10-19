# MED12–SUPT6H–H2A.Z Analysis

Analysis scripts for investigating the MED12–SUPT6H axis, H2A.Z-associated chromatin organization, transcriptional pause release, and gene expression.

[中文说明](README.zh-CN.md)

## Overview

This repository brings together scripts for APEX proteomics, ChIP-seq, PRO-seq, TT-seq, MNase-seq, and mouse liver RNA-seq. It includes differential analysis, genomic signal quantification, pathway enrichment, and integration of chromatin occupancy with transcriptional responses.

The repository currently provides a curated collection of existing scripts. Some require local path configuration and external input files; a complete end-to-end workflow has not been validated.

## Repository structure

```text
archive/
├── MS/                 # APEX proteomics
├── CHIP-SEQ/           # ChIP-seq annotation and regional signals
├── pro/                # PRO-seq pausing analysis
├── TT/                 # TT-seq analysis and metagene preparation
├── MNase/              # Nucleosome positional analysis (Python)
├── mouse/              # Mouse RNA-seq differential analysis
├── mRNAseq/            # H2A.Z and RNA-seq integration
└── med12_disease/      # External dataset comparisons
docs/                  # Input requirements, known issues, and validation
SOURCE_MANIFEST.tsv     # Original source locations and SHA-256 checksums
requirements-analysis.txt
```

Analysis scripts are organized under `archive/`. The MNase module is a Python conversion; other retained scripts preserve their source content. See the source manifest for provenance.

## Requirements

- Python and the packages listed in `requirements-analysis.txt`.
- R with module-specific packages, including DESeq2, tidyverse, clusterProfiler, and ChIPseeker where applicable.
- Linux or WSL for shell workflows and their external sequencing tools.
- Compatible genome annotations, input tables, and sequencing coverage files.
- Arial for plotting scripts that explicitly require it.

See [input requirements and execution notes](docs/INPUTS_AND_EXECUTION.md) for module-specific dependencies. The dependency list is not a tested environment lockfile.

## Getting started

Clone the repository and install Python dependencies in a dedicated environment:

```bash
git clone https://github.com/Berialim/MED12-SUPT6H-H2AZ-analysis.git
cd MED12-SUPT6H-H2AZ-analysis
python -m venv .venv
# Activate .venv using the command appropriate for your shell.
python -m pip install -r requirements-analysis.txt
```

Select the relevant analysis module, prepare its required inputs, and review file paths and parameters before execution. R packages and command-line tools must be installed separately for the modules that use them.

## Usage

The following scripts accept explicit input paths:

```bash
python archive/TT/plot_h2az_boxplot.py /path/to/results.txt \
  --density /path/to/density.bed -o /path/to/tt_h2az.pdf

python archive/pro/plot_pauseindex_h2az_boxplot.py /path/to/pauseindex.tsv \
  --density /path/to/density.bed \
  --ctrl-pause CTRL_pause --kd-pause KD_pause \
  -o /path/to/pro_h2az.pdf
```

The density input uses column 4 for gene identifiers and column 7 for H2A.Z signal by default. The PRO-seq plotting input must contain a `gene` column and the specified pause-index columns. These examples visualize prepared results; they do not run raw-read processing.

Other archived scripts may contain fixed paths or working-directory assumptions. Inspect them before use and run processing scripts on working copies of data.

## Data availability

The associated manuscript lists GEO accessions **GSE325384**, **GSE325382**, **GSE325378**, and **GSE325610**. Public availability and the mapping between accessions and assays have not been independently verified here.

Raw sequencing data, reference genomes, proteomics input tables, and gene-set collections are not included in this repository.

## Documentation and validation

- [Inputs and execution](docs/INPUTS_AND_EXECUTION.md)
- [Known issues and limitations](docs/KNOWN_ISSUES.md)
- [Validation report](docs/VALIDATION.md)
- [Source manifest](SOURCE_MANIFEST.tsv)

Checks cover source-file integrity and Python syntax. Scientific analyses have not been rerun. Known issues include strand handling in the archived pausing calculation and inconsistent grouping in the MNase Python script; review the linked documentation before interpreting outputs.

## Citation

Associated manuscript: *MED12 recruits SUPT6H to H2A.Z nucleosomes for pause release*.

Final publication details and a DOI will be added when available.

## License

No license has been assigned. A reuse license remains to be specified by the authors.

MNase usage: `python archive/MNase/mnase_stability.py /path/to/working_directory`.
