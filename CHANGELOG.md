# Changelog

All notable public changes to ProteoGenScope are documented in this file.

The project follows tool-level versioning: internal experimental pipeline revisions are not exposed as independent public software versions.

## [2.0] - 2026-09-25

### Added

- Public `proteoform-oriented` profile.
- Sample-specific nested DB1 / DB2 / DB3 search-space construction.
- EMBOSS `getorf` START-to-STOP candidate generation with a 390 nt minimum ORF size.
- Sense and reverse-complement ORF enumeration with explicit orientation provenance.
- TransDecoder as an additional coding-evidence branch.
- Exact full-length amino-acid sequence non-redundancy using SHA-256 sequence keys.
- Stable sample-specific PGS accessions.
- Proteoform registry, provenance, orientation, database membership and summary outputs.
- Sanitized samplesheet and metadata examples for public reproducibility.
- SLURM submission workflow for one independent job per biological sample.

### Validation

- Executed on three HNSCC RNA-seq samples in the CNPEM Marvin HPC environment.
- DB1 contained approximately 16.4k–17.8k non-redundant RNA-derived candidates per sample.
- DB1 candidates had a minimum length of 130 aa and a mean length of approximately 205 aa.

## Earlier development

Earlier pipeline revisions were internal development iterations and are intentionally not represented as separate public ProteoGenScope software versions.
