#!/usr/bin/env bash
set -euo pipefail

# Run in two stages so the profile matrices can be inspected before PDF creation:
#   run_ctrl_metagene.sh compute
#   run_ctrl_metagene.sh plot
ROOT=/mnt/d/work/pby
DT=/home/beria/.cache/pby-deeptools/venv
PLOT_CONFIG="$ROOT/TT/metagene_plot_config"
MATRIX_DIR="$ROOT/TT/.tmp_deeptools"

case "${1:-}" in
  compute)
    mkdir -p "$MATRIX_DIR"
    "$DT/bin/python" "$ROOT/TT/prepare_metagene_beds.py"

    "$DT/bin/computeMatrix" reference-point \
      --referencePoint TSS -b 100 -a 800 --binSize 10 \
      --missingDataAsZero --sortRegions keep -p 4 \
      -R "$ROOT/TT/med12/med12_PADJ0.05_UP_TSS.bed" \
         "$ROOT/TT/med12/med12_PADJ0.05_DOWN_TSS.bed" \
      -S "$ROOT/MNase/CTRL.bw" \
      -o "$MATRIX_DIR/med12_ctrl_tss_matrix.gz"

    "$DT/bin/computeMatrix" reference-point \
      --referencePoint TSS -b 100 -a 800 --binSize 10 \
      --missingDataAsZero --sortRegions keep -p 4 \
      -R "$ROOT/TT/spt6/spt6_PADJ0.05_UP_TSS.bed" \
         "$ROOT/TT/spt6/spt6_PADJ0.05_DOWN_TSS.bed" \
      -S "$ROOT/MNase/CTRL.bw" \
      -o "$MATRIX_DIR/spt6_ctrl_tss_matrix.gz"
    ;;
  plot)
    test -s "$MATRIX_DIR/med12_ctrl_tss_matrix.gz"
    test -s "$MATRIX_DIR/spt6_ctrl_tss_matrix.gz"
    MED12_UP_N=$(wc -l < "$ROOT/TT/med12/med12_PADJ0.05_UP_TSS.bed")
    MED12_DOWN_N=$(wc -l < "$ROOT/TT/med12/med12_PADJ0.05_DOWN_TSS.bed")
    SPT6_UP_N=$(wc -l < "$ROOT/TT/spt6/spt6_PADJ0.05_UP_TSS.bed")
    SPT6_DOWN_N=$(wc -l < "$ROOT/TT/spt6/spt6_PADJ0.05_DOWN_TSS.bed")

    FONTCONFIG_FILE="$PLOT_CONFIG/fonts.conf" \
    MPLCONFIGDIR="$PLOT_CONFIG/mplconfig" \
      "$DT/bin/plotProfile" \
        -m "$MATRIX_DIR/med12_ctrl_tss_matrix.gz" \
        -out "$ROOT/TT/med12/med12_PADJ0.05_UP_vs_DOWN_CTRL_TSS_metagene.pdf" \
        --plotType lines --averageType mean \
        --samplesLabel CTRL \
        --regionsLabel "UP (n=${MED12_UP_N})" "DOWN (n=${MED12_DOWN_N})" \
        --plotTitle "MED12 padj < 0.05: UP vs DOWN; TSS -100 to +800 bp" \
        --refPointLabel "TSS (0 bp)" \
        --yAxisLabel "Mean CTRL BigWig signal" \
        --plotWidth 18 --plotHeight 12 --colors "#0072B2" "#D55E00" \
        --plotFileFormat pdf

    FONTCONFIG_FILE="$PLOT_CONFIG/fonts.conf" \
    MPLCONFIGDIR="$PLOT_CONFIG/mplconfig" \
      "$DT/bin/plotProfile" \
        -m "$MATRIX_DIR/spt6_ctrl_tss_matrix.gz" \
        -out "$ROOT/TT/spt6/spt6_PADJ0.05_UP_vs_DOWN_CTRL_TSS_metagene.pdf" \
        --plotType lines --averageType mean \
        --samplesLabel CTRL \
        --regionsLabel "UP (n=${SPT6_UP_N})" "DOWN (n=${SPT6_DOWN_N})" \
        --plotTitle "SPT6 padj < 0.05: UP vs DOWN; TSS -100 to +800 bp" \
        --refPointLabel "TSS (0 bp)" \
        --yAxisLabel "Mean CTRL BigWig signal" \
        --plotWidth 18 --plotHeight 12 --colors "#0072B2" "#D55E00" \
        --plotFileFormat pdf
    rm -rf "$MATRIX_DIR"
    ;;
  *)
    echo "Usage: $0 {compute|plot}" >&2
    exit 2
    ;;
esac
