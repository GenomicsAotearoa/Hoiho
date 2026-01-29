#!/usr/bin/env bash

# End-to-end BCF processing pipeline.
# Steps:
#  1) Remove specific samples
#  2) Keep SNPs only and biallelic (m2/M2)
#  3) Filter QUAL >= 20
#  4) Annotate F_MISSING and filter F_MISSING <= 0.2
#  5) Split by populations (reusing split_populations.bash) [split-only mode]
#  6) Per-population MAC >= 2 filtering
#  7) At each step: write indexed BCF and bcftools stats; print SNP counts
#
# Usage:
#   ./process_bcf.bash [input.bcf] [pop_assignments.txt]
# Defaults:
#   input.bcf = merged.a9.unfiltered.bcf
#   pop_assignments.txt = pop_assignments.txt
#
# Filenames include a date tag. By default this script uses the requested date
# string for 8 Sept 2025: 2025-09-08. You can override via DATE_TAG env var.

set -euo pipefail

# Resolve bcftools command; fall back to pixi if needed.
if command -v bcftools >/dev/null 2>&1; then
  BCFTOOLS="bcftools"
else
  BCFTOOLS="pixi run bcftools"
fi

INPUT_BCF="${1:-merged.a9.unfiltered.bcf}"
POP_FILE="${2:-pop_assignments.txt}"
DATE_TAG="${DATE_TAG:-2025-09-08}"
THREADS="${THREADS:-0}"

if [[ ! -f "$INPUT_BCF" ]]; then
  echo "ERROR: Input BCF not found: $INPUT_BCF" >&2
  exit 1
fi
if [[ ! -f "$POP_FILE" ]]; then
  echo "ERROR: Population assignment file not found: $POP_FILE" >&2
  exit 1
fi

base_noext="${INPUT_BCF%.bcf}"

log() { printf "[%s] %s\n" "$(date +%H:%M:%S)" "$*"; }

run_index() {
  local bcf="$1"
  $BCFTOOLS index -f "$bcf"
}

run_stats() {
  local label="$1" bcf="$2" outf="$3"
  $BCFTOOLS stats "$bcf" > "$outf"
  local snps
  snps=$(awk -F"\t" '$1=="SN" && $3=="number of SNPs:" {print $4}' "$outf" | tail -n1)
  if [[ -z "$snps" ]]; then snps="NA"; fi
  printf "SNPs after %-22s: %s\n" "$label" "$snps"
}

# Stage 0: optionally report baseline stats
baseline_stats="${base_noext}.${DATE_TAG}.baseline.stats.txt"
run_stats "baseline" "$INPUT_BCF" "$baseline_stats"

# Stage 1: Remove specific samples
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT
cat > "$tmpdir/exclude_samples.txt" <<EOF
PP6
PP19
N23-25
YEP14
EOF

out1="${base_noext}.${DATE_TAG}.removed_samples.bcf"
log "Removing samples: PP6, PP19, N23-25, YEP14 -> $out1"
$BCFTOOLS view -S ^"$tmpdir/exclude_samples.txt" -Ob -o "$out1" --threads "$THREADS" "$INPUT_BCF"
run_index "$out1"
run_stats "remove_samples" "$out1" "${out1%.bcf}.stats.txt"

# Stage 2: SNPs only and biallelic
out2="${base_noext}.${DATE_TAG}.removed_samples.snps_biallelic.bcf"
log "Filtering to SNPs only and biallelic -> $out2"
$BCFTOOLS view -v snps -m2 -M2 -Ob -o "$out2" --threads "$THREADS" "$out1"
run_index "$out2"
run_stats "snps_biallelic" "$out2" "${out2%.bcf}.stats.txt"

# Stage 3: QUAL >= 20
out3="${base_noext}.${DATE_TAG}.removed_samples.snps_biallelic.qualge20.bcf"
log "Filtering QUAL >= 20 -> $out3"
$BCFTOOLS view -i 'QUAL>=20' -Ob -o "$out3" --threads "$THREADS" "$out2"
run_index "$out3"
run_stats "qual_ge_20" "$out3" "${out3%.bcf}.stats.txt"

# Stage 4: F_MISSING <= 0.2 (ensure tag present via +fill-tags)
out4="${base_noext}.${DATE_TAG}.removed_samples.snps_biallelic.qualge20.fmissing_le0.2.bcf"
log "Annotating F_MISSING and filtering <= 0.2 -> $out4"
$BCFTOOLS +fill-tags "$out3" -- -t F_MISSING | \
  $BCFTOOLS view -i 'F_MISSING<=0.2' -Ob -o "$out4" --threads "$THREADS"
run_index "$out4"
run_stats "f_missing_le_0.2" "$out4" "${out4%.bcf}.stats.txt"

# Stage 5: Split by populations (split-only, preserve DATE_TAG in outputs)
log "Splitting by populations using split_populations.bash (split-only)"
SPLIT_ONLY=1 OUT_TAG="split_${DATE_TAG}" THREADS="$THREADS" ./split_populations.bash "$out4" "$POP_FILE"

# The split script outputs: <base>.split_<DATE_TAG>.{northern,enderby,campbell}.bcf
declare -a POPS=(northern enderby campbell)
for pop in "${POPS[@]}"; do
  pop_bcf="${out4%.bcf}.split_${DATE_TAG}.${pop}.bcf"
  if [[ ! -f "$pop_bcf" ]]; then
    echo "WARNING: Expected split BCF not found: $pop_bcf" >&2
    continue
  fi
  run_index "$pop_bcf"  # ensure indexed (split script already indexes)
  run_stats "split_${pop}" "$pop_bcf" "${pop_bcf%.bcf}.stats.txt"

  # Stage 6: Apply MAC >= 2 per population
  pop_mac_bcf="${pop_bcf%.bcf}.macge2.bcf"
  log "Applying MAC >= 2 to ${pop} -> $pop_mac_bcf"
  $BCFTOOLS view -i 'MAC>=2' -Ob -o "$pop_mac_bcf" --threads "$THREADS" "$pop_bcf"
  run_index "$pop_mac_bcf"
  run_stats "mac_ge_2_${pop}" "$pop_mac_bcf" "${pop_mac_bcf%.bcf}.stats.txt"
done

log "All steps complete. Outputs are BCFs with .csi and .stats.txt at each stage."

