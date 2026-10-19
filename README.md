# MED12 SUPT6H H2A.Z analysis

Code collection supporting the manuscript **MED12 recruits SUPT6H to H2A.Z nucleosomes for pause release** (manuscript dated 17 October 2026).

This repository organizes existing local analysis scripts for APEX proteomics, ChIP-seq, PRO-seq, TT-seq, MNase-seq and mouse liver RNA-seq. It is a partial analysis-code archive, not a validated end-to-end reproduction of the manuscript. Figure associations are candidates established from code content and manuscript methods; final figure-generation provenance has not been established.

## Analysis modules

| Module | Manuscript relationship | Location |
| --- | --- | --- |
| APEX proteomics | Fig. 1 and Supplementary Fig. 1 candidates | `archive/MS/` |
| ChIP-seq processing and annotation | Fig. 2 and Supplementary Fig. 2 candidates | `archive/CHIP-SEQ/`, `archive/chipseq_se.sh`, `archive/spikein.sh` |
| PRO-seq pausing | Fig. 3b–c and Fig. 4i candidates | `archive/pro/` |
| TT-seq differential transcription | Fig. 3g–h and Fig. 4j candidates | `archive/TT/` |
| MNase positional stability | Supplementary Fig. 5c candidate | `notebooks/`, `archive/MNase/` |
| Metagene visualization | Supporting plotting utility | `archive/ratio/` |
| Mouse RNA-seq and H2A.Z integration | Fig. 6 and Supplementary Fig. 7 candidates | `archive/mouse/`, `archive/mRNAseq/` |
| External MED12 comparison | Supplementary Fig. 7 candidate | `archive/med12_disease/` |

## Start here

1. Read [known issues](docs/KNOWN_ISSUES.md) before executing archived code.
2. Consult [inputs and execution](docs/INPUTS_AND_EXECUTION.md) and restore required data/reference files.
3. Use `SOURCE_MANIFEST.tsv` for source locations and SHA-256 checksums. All files under `archive/` are byte-identical source copies.
4. Inspect [validation](docs/VALIDATION.md) for what has and has not been checked.

Two plotting scripts already accept explicit inputs (run from repository root after preparing data and dependencies):

```bash
python archive/TT/plot_h2az_boxplot.py /path/to/results.txt --density /path/to/density.bed -o /path/to/tt_h2az.pdf
python archive/pro/plot_pauseindex_h2az_boxplot.py /path/to/pauseindex.tsv --density /path/to/density.bed --ctrl-pause CTRL_pause --kd-pause KD_pause -o /path/to/pro_h2az.pdf
```

These commands are usage examples, not evidence of successful reproduction. Arial is required by these plotting scripts. Archived scripts retain original paths and assumptions. No script runs automatically.

## Data availability

The manuscript lists GEO accessions GSE325384, GSE325382, GSE325378 and GSE325610. Their public availability and assay-to-accession mapping were not independently checked. Sequencing files, manuscript text, proteomics tables, patient-level data and reference gene sets are not bundled.

## License and citation

No redistribution license is assigned by this archive. The authors must choose a license and confirm rights before public release. Cite the manuscript using its final publication metadata when available; no DOI is asserted here.
