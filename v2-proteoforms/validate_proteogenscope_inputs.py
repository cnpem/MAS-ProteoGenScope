#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path

SHEET_COLS = ["dataset_id", "sample_id", "r1", "r2"]
META_COLS = [
    "dataset_id", "patient_id", "sample_id", "original_sample_label",
    "sample_class", "tnm_label", "tnm_basis", "t_category", "n_category",
    "m_category", "nodal_status",
]


def read_tsv(path, required):
    with Path(path).open("r", encoding="utf-8", newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        if not r.fieldnames:
            raise SystemExit(f"ERROR: no header: {path}")
        missing = [x for x in required if x not in r.fieldnames]
        if missing:
            raise SystemExit(f"ERROR: {path} missing columns: {', '.join(missing)}")
        return list(r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--check-fastq", action="store_true")
    a = ap.parse_args()

    samples = read_tsv(a.samplesheet, SHEET_COLS)
    metadata = read_tsv(a.metadata, META_COLS)
    if not samples:
        raise SystemExit("ERROR: samplesheet contains no samples")

    meta_by_id = {}
    for row in metadata:
        sid = row["sample_id"]
        if sid in meta_by_id:
            raise SystemExit(f"ERROR: duplicated sample_id in metadata: {sid}")
        meta_by_id[sid] = row

    seen = set()
    for row in samples:
        sid = row["sample_id"]
        if sid in seen:
            raise SystemExit(f"ERROR: duplicated sample_id in samplesheet: {sid}")
        seen.add(sid)
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", sid):
            raise SystemExit(f"ERROR: invalid sample_id: {sid}")
        if sid not in meta_by_id:
            raise SystemExit(f"ERROR: sample_id absent from metadata: {sid}")
        if row["dataset_id"] != meta_by_id[sid]["dataset_id"]:
            raise SystemExit(
                f"ERROR: dataset mismatch for {sid}: samplesheet={row['dataset_id']} "
                f"metadata={meta_by_id[sid]['dataset_id']}"
            )
        if a.check_fastq:
            for col in ("r1", "r2"):
                p = Path(row[col])
                if not p.is_file() or p.stat().st_size == 0:
                    raise SystemExit(f"ERROR: missing/empty {col} for {sid}: {p}")

    print(f"OK: {len(samples)} samples validated")
    for row in samples:
        m = meta_by_id[row["sample_id"]]
        print(
            f"  {row['sample_id']} | patient={m['patient_id']} | "
            f"TNM={m['tnm_label']} | basis={m['tnm_basis']}"
        )


if __name__ == "__main__":
    main()
