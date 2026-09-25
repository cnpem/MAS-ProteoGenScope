#!/bin/bash
set -euo pipefail

PGS_ROOT="${PROTEOGENSCOPE_ROOT:?Set PROTEOGENSCOPE_ROOT}"
SAMPLESHEET="${SAMPLESHEET:-${PGS_ROOT}/samplesheet_proteogenscope.tsv}"
METADATA="${METADATA:-${PGS_ROOT}/metadata_proteogenscope.tsv}"
PIPELINE="${PGS_ROOT}/ProteogenScope_v2_proteoforms.sbatch"
VALIDATOR="${PGS_ROOT}/validate_proteogenscope_inputs.py"

CANONICAL="${UNIPROT_CANONICAL:?Set UNIPROT_CANONICAL}"
ISOFORMS="${UNIPROT_ISOFORMS:?Set UNIPROT_ISOFORMS}"
OUTROOT="${OUTROOT:?Set OUTROOT}"

# Proteoform profile convention.
ORF_MODE="start2stop"
ORF_MIN_NT=390
ORF_REVERSE="Y"

python3 "$VALIDATOR" --samplesheet "$SAMPLESHEET" --metadata "$METADATA" --check-fastq
for f in "$PIPELINE" "$CANONICAL" "$ISOFORMS"; do
  [[ -s "$f" ]] || { echo "ERROR: missing/empty file: $f" >&2; exit 1; }
done
mkdir -p "$OUTROOT/logs"

printf "%-18s %-12s %s\n" "SAMPLE_ID" "JOB_ID" "OUTDIR"
printf "%-18s %-12s %s\n" "------------------" "------------" "------"

while IFS=$'\t' read -r DATASET_ID SAMPLE_ID R1 R2; do
  [[ "$DATASET_ID" == "dataset_id" ]] && continue
  [[ -z "$DATASET_ID" ]] && continue
  SAMPLE_OUT="${OUTROOT}/${DATASET_ID}/${SAMPLE_ID}"
  mkdir -p "$SAMPLE_OUT"
  JOB_ID=$(sbatch --parsable \
    --job-name="PGS7PF_${SAMPLE_ID}" \
    --output="${OUTROOT}/logs/${SAMPLE_ID}_%j.out" \
    --error="${OUTROOT}/logs/${SAMPLE_ID}_%j.err" \
    "$PIPELINE" \
      --sample-id "$SAMPLE_ID" \
      --r1 "$R1" \
      --r2 "$R2" \
      --metadata "$METADATA" \
      --canonical "$CANONICAL" \
      --isoforms "$ISOFORMS" \
      --orf-mode "$ORF_MODE" \
      --orf-min-nt "$ORF_MIN_NT" \
      --orf-reverse "$ORF_REVERSE" \
      --outdir "$SAMPLE_OUT")
  printf "%-18s %-12s %s\n" "$SAMPLE_ID" "$JOB_ID" "$SAMPLE_OUT"
done < "$SAMPLESHEET"

echo
echo "ProteogenScope v2 proteoform jobs submitted."
echo "Profile: ${ORF_MODE}, min=${ORF_MIN_NT} nt, reverse=${ORF_REVERSE}"
echo "Monitor: squeue -u \"$USER\""
echo "Logs: $OUTROOT/logs"
