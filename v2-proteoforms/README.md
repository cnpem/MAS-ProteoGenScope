# ProteoGenScope v2 — Proteoform-oriented profile

This directory contains the current **proteoform-oriented implementation of ProteoGenScope v2**.

ProteoGenScope versioning refers to the tool as a whole. The term `proteoforms` here identifies the analytical profile implemented by this workflow and should not be interpreted as a separate software version.

## Purpose

The profile builds sample-specific proteogenomic FASTA search spaces from RNA-seq-supported transcripts while preserving the origin and orientation of each candidate protein sequence.

The workflow separates:

1. RNA-supported transcript reconstruction;
2. ORF candidate generation;
3. coding-evidence annotation;
4. exact amino-acid sequence non-redundancy;
5. incremental construction of DB1, DB2 and DB3;
6. provenance and membership reporting.

## ORF candidate convention

RNA-supported transcript sequences are processed with EMBOSS `getorf` using:

```text
-find 1        START -> STOP
-minsize 390   >=390 nt (~130 aa)
-reverse Y     enumerate sense and reverse-complement candidates
-table 0       standard nuclear genetic code
```

The current policy therefore generates START-to-STOP ORFs of at least 130 aa from both orientations.

`reverse=Y` is retained deliberately. Rather than discarding the distinction between orientations, ProteoGenScope propagates orientation into the provenance records so that sense and reverse-derived candidates can be evaluated separately downstream.

## Exact non-redundancy

Two candidates are collapsed only when their **complete amino-acid sequences are identical**.

SHA-256 hashes of the normalized amino-acid sequences are used as deterministic sequence identity keys.

No approximate sequence clustering is performed.

Consequently:

```text
100% identical full-length AA sequence -> one sequence entity
any AA substitution, indel or length difference -> distinct sequence entity
```

If multiple source ORFs produce the same sequence, only one physical sequence is written to the FASTA while all source records remain represented in `proteoform_provenance.tsv`.

## Stable ProteoGenScope accessions

Sample-specific sequence entities receive accessions using:

```text
<sample_id>_PGS<6-digit integer>
```

Example:

```text
SAMPLE01_PGS000143
```

The same sequence entity retains the same accession across DB1, DB2 and DB3 for a given sample.

## Orientation and provenance

For RNA-derived candidates, ProteoGenScope records available EMBOSS-derived information including:

- `transcript_id`;
- `orf_index`;
- `orientation` (`sense` or `reverse`);
- `emboss_start_nt`;
- `emboss_end_nt`;
- amino-acid sequence length.

After exact-AA collapse, a sequence can be classified as:

```text
sense
reverse
sense+reverse
```

The last category indicates that the same non-redundant amino-acid sequence has provenance from ORFs detected in both orientations.

## TransDecoder evidence branch

`TransDecoder.LongOrfs` and `TransDecoder.Predict` are executed as an additional coding-evidence layer.

In ProteoGenScope v2, TransDecoder does **not** replace EMBOSS `getorf` as the candidate generator and is not a hard filter for DB1.

Exact sequence agreement with TransDecoder can instead be recorded as supporting evidence for a candidate sequence.

## Database hierarchy

The profile builds three nested databases:

```text
DB1 = NR(RNA-seq proteoform candidates)
DB2 = NR(DB1 + UniProt canonical)
DB3 = NR(DB2 + UniProt isoforms)
```

Every union is followed by exact full-length amino-acid sequence non-redundancy.

This yields three progressively expanded search spaces while retaining a common ProteoGenScope sequence registry.

## Main outputs

```text
proteogenscope_databases/
├── <sample>_DB1_RNAseq_Proteoforms_NR.fasta
├── <sample>_DB2_RNAseq_Proteoforms_UniProtCanonical_NR.fasta
├── <sample>_DB3_RNAseq_Proteoforms_UniProtCanonical_Isoforms_NR.fasta
├── proteoform_registry.tsv
├── proteoform_provenance.tsv
├── orientation_summary.tsv
├── db_membership.tsv
├── duplicate_groups.tsv
├── database_summary.tsv
├── source_summary.tsv
└── database_manifest.tsv
```

### `proteoform_registry.tsv`

Master table of non-redundant amino-acid sequence entities. It links each PGS accession to its sequence hash, length, source classes, orientation summary, supporting evidence and DB membership.

### `proteoform_provenance.tsv`

One-to-many record of the original sequence sources contributing to each non-redundant PGS entity. RNA-derived entries retain transcript/ORF identifiers, orientation and EMBOSS coordinates where available.

### `orientation_summary.tsv`

Summarizes sense and reverse contributions before and after exact sequence collapse.

### `database_summary.tsv`

Reports database-level statistics including unique sequence count, newly introduced sequences, total amino-acid count and minimum/mean/maximum sequence length.

### Membership and provenance support files

`db_membership.tsv`, `duplicate_groups.tsv`, `source_summary.tsv` and `database_manifest.tsv` provide the supporting records needed to reconstruct how each FASTA search database was assembled.

## FASTA headers

ProteoGenScope-specific FASTA records retain a stable PGS accession and compact provenance information.

Example:

```text
>SAMPLE01_PGS000143 db=DB1 origin=RNAseq orientation=reverse
```

## Running on SLURM

The public repository intentionally excludes real sample metadata, internal HPC paths and generated databases.

Start from the sanitized templates under `examples/`, configure the required environment variables, and submit the workflow:

```bash
export PROTEOGENSCOPE_ROOT=/path/to/MAS-ProteoGenScope/v2-proteoforms
export PGS_ENV_BIN=/path/to/conda/env/bin
export GENOME_DIR=/path/to/STAR/index
export GTF_REF=/path/to/reference.gtf
export STRINGTIE_BIN=/path/to/stringtie
export FILTER_SJ=/path/to/filter_sj.py
export ANNOTATE_JUNC=/path/to/annotate_junctions.py
export FILTRAR_GTF_FPKM=/path/to/filtrar_gtf_fpkm.py
export FILTRAR_GTF_SALMON=/path/to/filtrar_gtf_salmon.py
export RELATORIO_FINAL=/path/to/gerar_relatorio_final.py
export UNIPROT_CANONICAL=/path/to/human_reviewed_canonical.fasta
export UNIPROT_ISOFORMS=/path/to/human_reviewed_isoforms_only.fasta
export OUTROOT=/path/to/results

bash submit_proteogenscope_v2_proteoforms.sh
```

The current HPC implementation also depends on project-specific helper scripts for transcript filtering, splice-junction annotation and reporting. These should be configured explicitly in the local environment before execution.

## Versioning policy

This directory belongs to **ProteoGenScope v2**.

Future changes in search-space strategy should preferably be represented as analytical profiles, configuration changes or documented feature additions. The software version should change when the framework itself changes materially, rather than whenever an internal experimental pipeline is refactored.
