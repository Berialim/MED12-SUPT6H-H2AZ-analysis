# Validation report

Updated 2026-10-19. The repository contains 21 retained source scripts and one derived MNase Python script. All retained original files match the SHA-256 checksums in SOURCE_MANIFEST.tsv. The MNase row identifies its notebook source and the checksum of the converted Python file.

All Python files passed syntax parsing. Scientific analyses were not rerun; R and shell files were not executed. No dependency environment or manuscript reproduction is certified by these checks.

## Cleanup

- Replaced the original and output-cleared MNase notebooks and duplicate code export with one Python script. Added an explicit working-directory argument and removed display-only expressions. Preserved the original scientific logic and documented grouping issue.
- Removed the legacy single-end ChIP-seq helper, incomplete spike-in alignment helper and malformed generic metagene helper. These did not provide the complete workflows described by the project.
- Original files outside this repository are unchanged. Earlier repository versions remain available in Git history.

## Python syntax checks

- `archive/med12_disease/comparison_med12_disease_vs_oe/plot_reference_pathways_nes.py`: PASS
- `archive/med12_disease/gsea_med12disease.py`: PASS
- `archive/MNase/mnase_stability.py`: PASS
- `archive/mRNAseq/analyze_three_variable_relationships.py`: PASS
- `archive/mRNAseq/h2az/h2azzscore+log2fc.py`: PASS
- `archive/mRNAseq/plot_score_20bin_barplots.py`: PASS
- `archive/mRNAseq/plot_score_relationships.py`: PASS
- `archive/MS/2023.09.24_MED12_APEX2/KO_MED12_CDK8_APEX2/gseapy_mediator_prerank_simple.py`: PASS
- `archive/MS/2023.09.24_MED12_APEX2/KO_MED12_CDK8_APEX2/volcano_mediator_down_simple.py`: PASS
- `archive/MS/2023.09.24_MED12_APEX2/KO_MED12_PRB1_APEX2/med12_prb1_simple_analysis.py`: PASS
- `archive/MS/2023.09.24_MED12_APEX2/MED_CORE_heatmap/med_core_cdk8_vs_rnapii_heatmap.py`: PASS
- `archive/pro/plot_pauseindex_h2az_boxplot.py`: PASS
- `archive/pro/spt6/pauseindex_proseq.py`: PASS
- `archive/TT/compairett.py`: PASS
- `archive/TT/plot_h2az_boxplot.py`: PASS
- `archive/TT/prepare_metagene_beds.py`: PASS
