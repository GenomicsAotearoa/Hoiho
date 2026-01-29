#!/usr/bin/env bash

# Filter out singletons (MAC >= 2) and/or split into three populations using bcftools.
# Usage:
#   ./split_populations.bash [input.bcf] [pop_assignments.txt]
# Defaults:
#   input.bcf = merged.a9.filtered.qual20_fmissing0.2.2alleles.snpsonly.bcf
#   pop_assignments.txt = pop_assignments.txt
#
# Behavior:
#  - By default, compute MAC (Minor Allele Count) and retain variants with MAC >= 2,
#    producing <input>.macge2.bcf and then split that BCF per population.
#  - If the environment variable SPLIT_ONLY=1 is set, the script will skip the MAC step
#    and only split the provided input file per population.
#  - You can customize the output tag inserted into output filenames by setting OUT_TAG.
#    Defaults: OUT_TAG="macge2" (normal mode) or OUT_TAG="split" (SPLIT_ONLY=1 mode).
#
# Outputs (default naming):
#  - Normal mode: <input>.macge2.{northern,enderby,campbell}.bcf (+ .csi)
#  - Split-only:  <input>.split.{northern,enderby,campbell}.bcf (+ .csi)

set -euo pipefail

# Resolve bcftools command; fall back to pixi if needed.
if command -v bcftools >/dev/null 2>&1; then
  BCFTOOLS="bcftools"
else
  BCFTOOLS="pixi run bcftools"
fi

INPUT_BCF="${1:-merged.a9.filtered.qual20_fmissing0.2.2alleles.snpsonly.bcf}"
POP_FILE="${2:-pop_assignments.txt}"
# Optional: number of compression threads for bcftools (0 = auto/disabled)
THREADS="${THREADS:-0}"
SPLIT_ONLY="${SPLIT_ONLY:-0}"

if [[ ! -f "$INPUT_BCF" ]]; then
  echo "ERROR: Input BCF not found: $INPUT_BCF" >&2
  exit 1
fi

if [[ ! -f "$POP_FILE" ]]; then
  echo "ERROR: Population assignment file not found: $POP_FILE" >&2
  exit 1
fi

BASE="${INPUT_BCF%.bcf}"

if [[ "$SPLIT_ONLY" != "1" ]]; then
  OUT_TAG="${OUT_TAG:-macge2}"
  FILTERED_BCF="${BASE}.${OUT_TAG}.bcf"
  echo "[1/3] Filtering to MAC>=2 via expression -> $FILTERED_BCF"
  # Use an explicit expression. MAC is computed on-the-fly if INFO/AC/AN absent.
  $BCFTOOLS view -i 'MAC>=2' -Ob -o "$FILTERED_BCF" --threads "$THREADS" "$INPUT_BCF"
  echo "[2/3] Indexing $FILTERED_BCF"
  $BCFTOOLS index -f "$FILTERED_BCF"
else
  OUT_TAG="${OUT_TAG:-split}"
  FILTERED_BCF="$INPUT_BCF"  # No MAC filtering in split-only mode
  echo "[1/2] Split-only mode: skipping MAC filter"
fi

# Prepare temporary directory for sample lists.
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

# Extract sample lists for the three populations.
declare -a POPS=("Northern" "Enderby" "Campbell")
for pop in "${POPS[@]}"; do
  sampfile="$tmpdir/${pop}.samples.txt"
  # Expecting two columns: ID and Population label; skip header line.
  awk -v target_pop="$pop" 'NR>1 && tolower($2)==tolower(target_pop) {print $1}' "$POP_FILE" > "$sampfile"
  if [[ ! -s "$sampfile" ]]; then
    echo "WARNING: No samples found for population: $pop" >&2
  fi
done

echo "[*/ ] Splitting by population"
for pop in "${POPS[@]}"; do
  sampfile="$tmpdir/${pop}.samples.txt"
  outbcf="${BASE}.${OUT_TAG}.${pop,,}.bcf"
  echo "  - $pop -> $outbcf"
  $BCFTOOLS view --force-samples -S "$sampfile" -Ob -o "$outbcf" --threads "$THREADS" "$FILTERED_BCF"
  $BCFTOOLS index -f "$outbcf"
done
echo "[Done] Outputs:"
if [[ "$SPLIT_ONLY" != "1" ]]; then
  echo "  - $FILTERED_BCF (+ .csi)"
fi
for pop in "${POPS[@]}"; do
  echo "  - ${BASE}.${OUT_TAG}.${pop,,}.bcf (+ .csi)"
done
