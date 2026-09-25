# ProteoGenScope v2 conventions

## Software version

```text
ProteoGenScope v2
```

The current implementation uses the analytical profile:

```text
proteoform-oriented
```

Internal development revisions are not exposed as separate public tool versions.

## Candidate ORF policy

- EMBOSS `getorf -find 1`: START-to-STOP ORFs;
- `-minsize 390`: minimum ORF size of 390 nt (~130 aa);
- `-reverse Y`: evaluate transcript sense and reverse complement;
- `-table 0`: standard nuclear genetic code.

## Sequence identity

A non-redundant ProteoGenScope sequence entity is defined by **exact full-length amino-acid sequence identity**.

```text
same complete AA sequence -> same sequence entity
any AA substitution / indel / length difference -> distinct sequence entity
```

SHA-256 hashes of normalized AA sequences are used as deterministic internal keys.

No similarity-based clustering threshold is used.

## Stable PGS accession

```text
<sample_id>_PGS<6-digit integer>
```

Example using a non-clinical identifier:

```text
SAMPLE01_PGS000143
```

Within one sample, the same sequence entity retains the same PGS accession across DB1, DB2 and DB3.

## Orientation

RNA-derived candidate provenance is classified as:

- `sense`;
- `reverse`;
- `sense+reverse` after exact-AA collapse when the same sequence has provenance from both orientations.

Orientation is evidence/provenance metadata and does not alter the sequence identity key.

## Database hierarchy

```text
DB1 = NR(RNA-seq-derived proteoform candidates)
DB2 = NR(DB1 + UniProt canonical)
DB3 = NR(DB2 + UniProt isoforms)
```

Every union is followed by exact full-length AA non-redundancy.

## TransDecoder

TransDecoder is treated as an additional coding-evidence source. It is not the sole candidate generator and is not currently a mandatory filter for DB1 admission.

## Provenance principle

Exact sequence collapse affects physical FASTA representation, not evidence retention.

If several source ORFs or reference records generate the same sequence:

```text
one AA sequence entity
many provenance records
```

This distinction is represented primarily by `proteoform_registry.tsv` and `proteoform_provenance.tsv`.
