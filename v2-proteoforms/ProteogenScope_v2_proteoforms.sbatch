#!/bin/bash
#SBATCH --job-name=ProteogenScope_PF
#SBATCH --partition=long-cpu
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=18
#SBATCH --mem=80G

set -euo pipefail

# ProteogenScope v2 — proteoform-oriented profile
# One sbatch submission = one biological sample = DB1/DB2/DB3.
# Default candidate policy: START->STOP, >=390 nt (130 aa), sense + reverse complement.

: "${PGS_ENV_BIN:?Set PGS_ENV_BIN to the bin directory of the ProteoGenScope environment}"
export PATH="${PGS_ENV_BIN}:$PATH"

SAMPLE_ID=""; R1=""; R2=""; METADATA=""; UNIPROT_CANONICAL=""; UNIPROT_ISOFORMS=""; OUTDIR=""; PREFIX=""
ORF_MODE="${ORF_MODE:-start2stop}"
ORF_MIN_NT="${ORF_MIN_NT:-390}"
ORF_REVERSE="${ORF_REVERSE:-Y}"

usage() {
  cat <<'USAGE'
Usage:
  sbatch ProteogenScope_v2_proteoforms.sbatch \
    --sample-id SAMPLE01 \
    --r1 sample_R1.fastq.gz \
    --r2 sample_R2.fastq.gz \
    --metadata metadata_proteogenscope_v0.2.tsv \
    --canonical human_reviewed_canonical.fasta \
    --isoforms human_reviewed_isoforms_only.fasta \
    --outdir /path/to/results/SAMPLE01

Proteoform defaults:
  --orf-mode start2stop
  --orf-min-nt 390
  --orf-reverse Y
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --r1) R1="$2"; shift 2 ;;
    --r2) R2="$2"; shift 2 ;;
    --metadata) METADATA="$2"; shift 2 ;;
    --canonical) UNIPROT_CANONICAL="$2"; shift 2 ;;
    --isoforms) UNIPROT_ISOFORMS="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --orf-mode) ORF_MODE="$2"; shift 2 ;;
    --orf-min-nt) ORF_MIN_NT="$2"; shift 2 ;;
    --orf-reverse) ORF_REVERSE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for v in SAMPLE_ID R1 R2 METADATA UNIPROT_CANONICAL UNIPROT_ISOFORMS OUTDIR; do
  [[ -n "${!v}" ]] || { echo "ERROR: missing required argument: $v" >&2; exit 2; }
done
[[ -n "$PREFIX" ]] || PREFIX="$SAMPLE_ID"

case "$ORF_MODE" in
  stop2stop) GETORF_FIND=0 ;;
  start2stop) GETORF_FIND=1 ;;
  *) echo "ERROR: --orf-mode must be stop2stop or start2stop" >&2; exit 2 ;;
esac
[[ "$ORF_MIN_NT" =~ ^[0-9]+$ ]] || { echo "ERROR: --orf-min-nt must be an integer" >&2; exit 2; }
ORF_REVERSE="${ORF_REVERSE^^}"
[[ "$ORF_REVERSE" == "Y" || "$ORF_REVERSE" == "N" ]] || { echo "ERROR: --orf-reverse must be Y or N" >&2; exit 2; }

THREADS="${SLURM_CPUS_PER_TASK:-18}"
: "${GENOME_DIR:?Set GENOME_DIR to the STAR genome-index directory}"
GENOME_FA="${GENOME_FA:-${GENOME_DIR}/genome.fa}"
GTF_REF="${GTF_REF:-${GENOME_DIR}/genes.gtf}"
STRINGTIE="${STRINGTIE_BIN:-stringtie}"
: "${FILTER_SJ:?Set FILTER_SJ to filter_sj.py}"
: "${ANNOTATE_JUNC:?Set ANNOTATE_JUNC to annotate_junctions.py}"
: "${FILTRAR_GTF_FPKM:?Set FILTRAR_GTF_FPKM to filtrar_gtf_fpkm.py}"
: "${FILTRAR_GTF_SALMON:?Set FILTRAR_GTF_SALMON to filtrar_gtf_salmon.py}"
: "${RELATORIO_FINAL:?Set RELATORIO_FINAL to gerar_relatorio_final.py}"
: "${PROTEOGENSCOPE_ROOT:?Set PROTEOGENSCOPE_ROOT to the v2-proteoforms repository directory}"
DB_BUILDER="${PROTEOGENSCOPE_ROOT}/build_proteogenscope_proteoforms_v2.py"
DB_DIR="${OUTDIR}/proteogenscope_databases"

mkdir -p "$OUTDIR" "$DB_DIR"
for f in "$R1" "$R2" "$GENOME_FA" "$GTF_REF" "$UNIPROT_CANONICAL" "$UNIPROT_ISOFORMS" "$METADATA" "$DB_BUILDER" "$FILTER_SJ" "$ANNOTATE_JUNC" "$FILTRAR_GTF_FPKM" "$FILTRAR_GTF_SALMON" "$RELATORIO_FINAL"; do
  [[ -s "$f" ]] || { echo "ERROR: missing/empty file: $f" >&2; exit 1; }
done
[[ -x "$STRINGTIE" ]] || { echo "ERROR: StringTie not executable: $STRINGTIE" >&2; exit 1; }
for cmd in STAR gffread gffcompare salmon getorf TransDecoder.LongOrfs TransDecoder.Predict python3; do
  command -v "$cmd" >/dev/null || { echo "ERROR: command not found: $cmd" >&2; exit 1; }
done

META_ROW=$(awk -F'\t' -v sid="$SAMPLE_ID" 'NR==1 {next} $3==sid {print; n++} END {if(n!=1) exit 42}' "$METADATA") || {
  echo "ERROR: sample_id must occur exactly once in metadata: $SAMPLE_ID" >&2; exit 1;
}

printf '%s\n' \
  "============================================================" \
  "ProteogenScope v2 — PROTEOFORMS" \
  "SLURM_JOB_ID : ${SLURM_JOB_ID:-NA}" \
  "Node         : $(hostname)" \
  "Sample       : $SAMPLE_ID" \
  "Output       : $OUTDIR" \
  "Threads      : $THREADS" \
  "ORF mode     : $ORF_MODE (getorf -find $GETORF_FIND)" \
  "ORF min nt   : $ORF_MIN_NT" \
  "ORF reverse  : $ORF_REVERSE" \
  "Metadata row : $META_ROW" \
  "============================================================"

cd "$OUTDIR"

echo "[1/14] STAR alignment"
STAR --runThreadN "$THREADS" --genomeDir "$GENOME_DIR" --readFilesIn "$R1" "$R2" --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "${PREFIX}_"

echo "[2/14] StringTie transcript assembly"
"$STRINGTIE" "${PREFIX}_Aligned.sortedByCoord.out.bam" -o "${PREFIX}_novo.gtf" -p "$THREADS" -l "${PREFIX}_novo"
"$STRINGTIE" "${PREFIX}_Aligned.sortedByCoord.out.bam" -G "$GTF_REF" -o "${PREFIX}.gtf" -p "$THREADS" -l "$PREFIX"

echo "[3/14] gffcompare"
gffcompare -r "${PREFIX}.gtf" -o "${PREFIX}_merged" "${PREFIX}_novo.gtf"

echo "[4/14] StringTie quantification"
"$STRINGTIE" "${PREFIX}_Aligned.sortedByCoord.out.bam" -e -B -G "${PREFIX}_merged.annotated.gtf" -o "${PREFIX}_quant.gtf" -p "$THREADS" -l "$PREFIX"

echo "[5/14] FPKM filter"
python3 "$FILTRAR_GTF_FPKM" "${PREFIX}_quant.gtf" "${PREFIX}_fpkm_filtered.gtf" --limiar_fpkm 1.0

echo "[6/14] Salmon quantification"
gffread "${PREFIX}_fpkm_filtered.gtf" -g "$GENOME_FA" -w "${PREFIX}_transcripts.fa"
salmon index -t "${PREFIX}_transcripts.fa" -i "${PREFIX}_salmon_index"
salmon quant -i "${PREFIX}_salmon_index" -l A --threads "$THREADS" -1 "$R1" -2 "$R2" -o "${PREFIX}_salmon_quant"
awk 'NR>1 && $4 > 1 {print $1}' "${PREFIX}_salmon_quant/quant.sf" > "${PREFIX}_expressed_transcripts.txt"

echo "[7/14] Final GTF curation"
python3 "$FILTRAR_GTF_SALMON" "${PREFIX}_fpkm_filtered.gtf" "${PREFIX}_expressed_transcripts.txt" "${PREFIX}_final_filtered.gtf"

echo "[8/14] Final RNA-supported transcripts"
gffread "${PREFIX}_final_filtered.gtf" -g "$GENOME_FA" -w "${PREFIX}_final_transcripts.fa"

echo "[9/14] Proteoform ORF candidate generation"
getorf \
  -sequence "${PREFIX}_final_transcripts.fa" \
  -outseq "${PREFIX}_proteoform_ORF_candidates.raw.fasta" \
  -find "$GETORF_FIND" \
  -minsize "$ORF_MIN_NT" \
  -table 0 \
  -reverse "$ORF_REVERSE" \
  -auto

echo "[10/14] TransDecoder evidence branch"
TransDecoder.LongOrfs -t "${PREFIX}_final_transcripts.fa"
TransDecoder.Predict -t "${PREFIX}_final_transcripts.fa"
TD_PEP="${PREFIX}_final_transcripts.fa.transdecoder.pep"

echo "[11/14] Build sample-specific proteoform DB1 / DB2 / DB3"
python3 "$DB_BUILDER" \
  --rnaseq "${PREFIX}_proteoform_ORF_candidates.raw.fasta" \
  --canonical "$UNIPROT_CANONICAL" \
  --isoforms "$UNIPROT_ISOFORMS" \
  --metadata "$METADATA" \
  --sample-id "$SAMPLE_ID" \
  --transdecoder-pep "$TD_PEP" \
  --orf-mode "$ORF_MODE" \
  --orf-min-nt "$ORF_MIN_NT" \
  --orf-reverse "$ORF_REVERSE" \
  --outdir "$DB_DIR"

echo "[12/14] Database + orientation summaries"
if command -v column >/dev/null 2>&1; then
  column -t -s $'\t' "$DB_DIR/database_summary.tsv"
  echo
  column -t -s $'\t' "$DB_DIR/orientation_summary.tsv"
else
  cat "$DB_DIR/database_summary.tsv"; echo; cat "$DB_DIR/orientation_summary.tsv"
fi

echo "[13/14] Junction analysis"
python3 "$FILTER_SJ" "${PREFIX}_SJ.out.tab" "${PREFIX}_filtered_junctions.csv" --min_unique_reads 10 --allowed_motifs 1 2
python3 "$ANNOTATE_JUNC" "${PREFIX}_filtered_junctions.csv" "${PREFIX}_junctions_annotated.csv"

echo "[14/14] Legacy HTML report"
python3 "$RELATORIO_FINAL" --quant "${PREFIX}_salmon_quant/quant.sf" --junctions "${PREFIX}_junctions_annotated.csv" --orfs "$TD_PEP" --output "${PREFIX}_relatorio_pipeline.html"

echo "DB1: $DB_DIR/${SAMPLE_ID}_DB1_RNAseq_Proteoforms_NR.fasta"
echo "DB2: $DB_DIR/${SAMPLE_ID}_DB2_RNAseq_Proteoforms_UniProtCanonical_NR.fasta"
echo "DB3: $DB_DIR/${SAMPLE_ID}_DB3_RNAseq_Proteoforms_UniProtCanonical_Isoforms_NR.fasta"
echo "DONE: $(date -Is)"
