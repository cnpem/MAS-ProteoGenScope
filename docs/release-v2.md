# ProteoGenScope v2 release notes

ProteoGenScope v2 formalizes the current **proteoform-oriented** search-space construction profile.

## Candidate generation

```text
EMBOSS getorf
-find 1
-minsize 390
-reverse Y
-table 0
```

This enumerates START-to-STOP ORFs of at least 390 nt (~130 aa) from both transcript orientations.

## Search-space hierarchy

```text
DB1 = NR(RNA-seq-derived proteoform candidates)
DB2 = NR(DB1 + UniProt canonical)
DB3 = NR(DB2 + UniProt isoforms)
```

Non-redundancy is defined by exact full-length amino-acid identity. SHA-256 sequence hashes are used as deterministic internal keys.

## Traceability

v2 records sequence identity separately from provenance and database membership. Sense/reverse orientation is retained for RNA-derived candidates, and TransDecoder is used as supporting coding evidence rather than as a mandatory DB1 filter.

## Validation snapshot

Across three HNSCC validation samples, DB1 contained approximately 16.4k–17.8k non-redundant RNA-derived candidates. Mean DB1 sequence length was ~205 aa, with a 130 aa minimum by design.
