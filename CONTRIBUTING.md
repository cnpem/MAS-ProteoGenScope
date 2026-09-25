# Contributing to ProteoGenScope

ProteoGenScope is under active development. Contributions that improve reproducibility, traceability, search-space control, portability or documentation are welcome.

## Before contributing

Please keep the following design principles intact unless the proposed change explicitly revises them:

1. candidate generation and search-space admission should remain conceptually separable;
2. exact sequence collapse must not erase provenance;
3. database membership should be reconstructable from tabular outputs;
4. public examples must not contain clinical identifiers, internal credentials or institution-specific private paths;
5. analytical assumptions should be explicit and configurable whenever possible.

## Suggested workflow

```bash
git clone https://github.com/cnpem/MAS-ProteoGenScope.git
cd MAS-ProteoGenScope
git checkout -b feature/short-description
```

Make focused changes and use descriptive commits. Before opening a pull request, verify shell and Python syntax where applicable:

```bash
bash -n v2-proteoforms/*.sh v2-proteoforms/*.sbatch
python -m py_compile v2-proteoforms/*.py
```

## Pull requests

A pull request should describe:

- the problem being addressed;
- the analytical or software behavior that changes;
- expected effects on DB1/DB2/DB3 construction;
- whether output schemas change;
- how the change was tested.

Changes that alter sequence identity, ORF admission, provenance or database membership should include a small reproducible example whenever possible.

## Data and privacy

Do not commit raw human sequencing data, patient-level metadata, generated sample databases, access tokens, credentials, internal-only URLs or private HPC paths.
