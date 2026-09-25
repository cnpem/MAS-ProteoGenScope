#!/usr/bin/env python3
# ProteogenScope database builder v2.0 — proteoform-oriented profile
from __future__ import annotations

import argparse
import csv
import hashlib
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Set

PROTEOGENSCOPE_VERSION = "2.0-proteoforms"
ENTITY_TYPE = "proteoform_candidate"

SOURCES = ("RNAseq", "UniProtCanonical", "UniProtIsoform")
REQUIRED_METADATA_COLUMNS = (
    "dataset_id", "patient_id", "sample_id", "original_sample_label",
    "sample_class", "tnm_label", "tnm_basis", "t_category",
    "n_category", "m_category", "nodal_status",
)
DB_ALLOWED_SOURCES = {
    "DB1": {"RNAseq"},
    "DB2": {"RNAseq", "UniProtCanonical"},
    "DB3": {"RNAseq", "UniProtCanonical", "UniProtIsoform"},
}
DB_DEFINITIONS = {
    "DB1": "NR(RNAseq proteoform candidates)",
    "DB2": "NR(RNAseq proteoform candidates + UniProtCanonical)",
    "DB3": "NR(RNAseq proteoform candidates + UniProtCanonical + UniProtIsoform)",
}


@dataclass(frozen=True)
class FastaRecord:
    source: str
    original_id: str
    original_header: str
    sequence: str
    orientation: str = "NA"
    transcript_id: str = "NA"
    orf_index: str = "NA"
    emboss_start_nt: str = "NA"
    emboss_end_nt: str = "NA"


def parse_getorf_header(header: str, source: str) -> tuple[str, str, str, str, str]:
    if source != "RNAseq":
        return "NA", "NA", "NA", "NA", "NA"

    upper = header.upper()
    orientation = "reverse" if "(REVERSE SENSE)" in upper else "sense"

    m_coord = re.search(r"\[\s*(\d+)\s*-\s*(\d+)\s*\]", header)
    if m_coord:
        start_nt, end_nt = m_coord.group(1), m_coord.group(2)
    else:
        start_nt, end_nt = "NA", "NA"

    original_id = header.split()[0]
    m_orf = re.match(r"^(.*)_([0-9]+)$", original_id)
    if m_orf:
        transcript_id, orf_index = m_orf.group(1), m_orf.group(2)
    else:
        transcript_id, orf_index = original_id, "NA"

    return orientation, transcript_id, orf_index, start_nt, end_nt


def read_fasta(path: Path, source: str) -> Iterator[FastaRecord]:
    header = None
    chunks: List[str] = []

    def emit() -> FastaRecord | None:
        nonlocal header, chunks
        if header is None:
            return None
        seq = "".join(chunks).replace(" ", "").replace("\t", "").upper().rstrip("*")
        if not seq:
            return None
        if "*" in seq:
            raise ValueError(f"Internal stop '*' in {source}: {header}")
        if not re.fullmatch(r"[A-Z]+", seq):
            raise ValueError(f"Invalid amino-acid characters in {source}: {header}")
        orientation, transcript_id, orf_index, start_nt, end_nt = parse_getorf_header(header, source)
        return FastaRecord(
            source=source,
            original_id=header.split()[0],
            original_header=header,
            sequence=seq,
            orientation=orientation,
            transcript_id=transcript_id,
            orf_index=orf_index,
            emboss_start_nt=start_nt,
            emboss_end_nt=end_nt,
        )

    with path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                rec = emit()
                if rec:
                    yield rec
                header = line[1:].strip()
                chunks = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                chunks.append(line)
        rec = emit()
        if rec:
            yield rec


def seq_hash(seq: str) -> str:
    return hashlib.sha256(seq.encode("ascii")).hexdigest()


def file_hash(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def validate_sample_id(sample_id: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample_id):
        raise SystemExit("Invalid --sample-id. Use only letters, numbers, underscore, dot and hyphen.")
    return sample_id


def load_metadata(path: Path, sample_id: str) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"Metadata has no header: {path}")
        missing = [c for c in REQUIRED_METADATA_COLUMNS if c not in reader.fieldnames]
        if missing:
            raise SystemExit(f"Metadata missing required columns: {', '.join(missing)}")
        matches = [row for row in reader if row.get("sample_id") == sample_id]
    if len(matches) != 1:
        raise SystemExit(f"sample_id={sample_id!r} must occur exactly once in metadata; found {len(matches)}")
    return {k: matches[0].get(k, "") for k in REQUIRED_METADATA_COLUMNS}


def introduced_in_db(h: str, source_hashes: Dict[str, Set[str]]) -> str:
    if h in source_hashes["RNAseq"]:
        return "DB1"
    if h in source_hashes["UniProtCanonical"]:
        return "DB2"
    if h in source_hashes["UniProtIsoform"]:
        return "DB3"
    raise RuntimeError(f"Sequence hash without source: {h}")


def origin_label(h: str, source_hashes: Dict[str, Set[str]], db: str) -> str:
    allowed = DB_ALLOWED_SOURCES[db]
    labels = [s for s in SOURCES if s in allowed and h in source_hashes[s]]
    return "+".join(labels)


def rnaseq_orientation_label(records: list[FastaRecord]) -> str:
    values = {r.orientation for r in records if r.source == "RNAseq" and r.orientation != "NA"}
    if not values:
        return "NA"
    ordered = [x for x in ("sense", "reverse") if x in values]
    return "+".join(ordered)


def write_fasta(
    path: Path,
    db: str,
    ordered_hashes: Iterable[str],
    members: Set[str],
    seqs: Dict[str, str],
    pgs_ids: Dict[str, str],
    source_hashes: Dict[str, Set[str]],
    provenance: dict[str, list[FastaRecord]],
) -> None:
    with path.open("w", encoding="utf-8") as out:
        for h in ordered_hashes:
            if h not in members:
                continue
            fields = [f"db={db}", f"origin={origin_label(h, source_hashes, db)}"]
            if h in source_hashes["RNAseq"]:
                fields.append(f"orientation={rnaseq_orientation_label(provenance[h])}")
            out.write(f">{pgs_ids[h]} {' '.join(fields)}\n")
            seq = seqs[h]
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Build sample-specific ProteogenScope v2 proteoform DB1/DB2/DB3 with exact-AA non-redundancy"
    )
    ap.add_argument("--rnaseq", required=True, type=Path)
    ap.add_argument("--canonical", required=True, type=Path)
    ap.add_argument("--isoforms", required=True, type=Path)
    ap.add_argument("--metadata", required=True, type=Path)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--transdecoder-pep", type=Path, default=None)
    ap.add_argument("--orf-mode", default="start2stop", choices=["stop2stop", "start2stop"])
    ap.add_argument("--orf-min-nt", type=int, default=390)
    ap.add_argument("--orf-reverse", default="Y", choices=["Y", "N", "y", "n"])
    args = ap.parse_args()

    sample_id = validate_sample_id(args.sample_id)
    for p in (args.rnaseq, args.canonical, args.isoforms, args.metadata):
        if not p.is_file() or p.stat().st_size == 0:
            raise SystemExit(f"Input file not found or empty: {p}")

    meta = load_metadata(args.metadata, sample_id)
    args.outdir.mkdir(parents=True, exist_ok=True)

    sources = [
        ("RNAseq", args.rnaseq),
        ("UniProtCanonical", args.canonical),
        ("UniProtIsoform", args.isoforms),
    ]

    seqs: Dict[str, str] = {}
    provenance: dict[str, list[FastaRecord]] = defaultdict(list)
    source_hashes: Dict[str, Set[str]] = defaultdict(set)
    raw_counts: dict[str, int] = defaultdict(int)

    for source, path in sources:
        for rec in read_fasta(path, source):
            raw_counts[source] += 1
            h = seq_hash(rec.sequence)
            if h in seqs and seqs[h] != rec.sequence:
                raise RuntimeError("SHA-256 collision detected")
            seqs[h] = rec.sequence
            source_hashes[source].add(h)
            provenance[h].append(rec)

    transdecoder_hashes: Set[str] = set()
    if args.transdecoder_pep and args.transdecoder_pep.is_file() and args.transdecoder_pep.stat().st_size:
        for rec in read_fasta(args.transdecoder_pep, "TransDecoder"):
            transdecoder_hashes.add(seq_hash(rec.sequence))

    db1 = set(source_hashes["RNAseq"])
    db2 = db1 | set(source_hashes["UniProtCanonical"])
    db3 = db2 | set(source_hashes["UniProtIsoform"])
    db_members = {"DB1": db1, "DB2": db2, "DB3": db3}

    rank = {"DB1": 1, "DB2": 2, "DB3": 3}
    ordered_hashes = sorted(seqs, key=lambda h: (rank[introduced_in_db(h, source_hashes)], h))
    width = max(6, len(str(len(ordered_hashes))))
    pgs_ids = {h: f"{sample_id}_PGS{i:0{width}d}" for i, h in enumerate(ordered_hashes, start=1)}

    fasta_paths = {
        "DB1": args.outdir / f"{sample_id}_DB1_RNAseq_Proteoforms_NR.fasta",
        "DB2": args.outdir / f"{sample_id}_DB2_RNAseq_Proteoforms_UniProtCanonical_NR.fasta",
        "DB3": args.outdir / f"{sample_id}_DB3_RNAseq_Proteoforms_UniProtCanonical_Isoforms_NR.fasta",
    }
    for db in ("DB1", "DB2", "DB3"):
        write_fasta(fasta_paths[db], db, ordered_hashes, db_members[db], seqs, pgs_ids, source_hashes, provenance)

    with (args.outdir / "sample_metadata.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(REQUIRED_METADATA_COLUMNS), delimiter="\t")
        w.writeheader(); w.writerow(meta)

    registry_path = args.outdir / "proteoform_registry.tsv"
    with registry_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "sample_id", "pgs_id", "entity_type", "sha256", "length_aa", "introduced_in_db",
            "from_rnaseq", "rnaseq_orientation", "n_rnaseq_sense", "n_rnaseq_reverse",
            "transdecoder_exact_support", "from_uniprot_canonical", "from_uniprot_isoform",
            "DB1", "DB2", "DB3", "n_provenance_records", "sequence",
        ])
        for h in ordered_hashes:
            recs = provenance[h]
            n_sense = sum(r.source == "RNAseq" and r.orientation == "sense" for r in recs)
            n_reverse = sum(r.source == "RNAseq" and r.orientation == "reverse" for r in recs)
            w.writerow([
                sample_id, pgs_ids[h], ENTITY_TYPE, h, len(seqs[h]), introduced_in_db(h, source_hashes),
                int(h in source_hashes["RNAseq"]), rnaseq_orientation_label(recs), n_sense, n_reverse,
                int(h in transdecoder_hashes), int(h in source_hashes["UniProtCanonical"]),
                int(h in source_hashes["UniProtIsoform"]), int(h in db1), int(h in db2), int(h in db3),
                len(recs), seqs[h],
            ])

    with (args.outdir / "proteoform_provenance.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "sample_id", "pgs_id", "source", "original_id", "transcript_id", "orf_index",
            "orientation", "emboss_start_nt", "emboss_end_nt", "aa_length", "original_header",
        ])
        for h in ordered_hashes:
            for rec in provenance[h]:
                w.writerow([
                    sample_id, pgs_ids[h], rec.source, rec.original_id, rec.transcript_id, rec.orf_index,
                    rec.orientation, rec.emboss_start_nt, rec.emboss_end_nt, len(rec.sequence), rec.original_header,
                ])

    with (args.outdir / "db_membership.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "pgs_id", "introduced_in_db", "DB1", "DB2", "DB3"])
        for h in ordered_hashes:
            w.writerow([sample_id, pgs_ids[h], introduced_in_db(h, source_hashes), int(h in db1), int(h in db2), int(h in db3)])

    with (args.outdir / "duplicate_groups.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "pgs_id", "sha256", "n_records", "sources", "orientations", "original_ids"])
        for h in ordered_hashes:
            recs = provenance[h]
            if len(recs) > 1:
                orientations = sorted({r.orientation for r in recs if r.orientation != "NA"})
                w.writerow([
                    sample_id, pgs_ids[h], h, len(recs), ";".join(sorted({r.source for r in recs})),
                    ";".join(orientations) if orientations else "NA", ";".join(r.original_id for r in recs),
                ])

    introduced_counts = {db: sum(1 for h in ordered_hashes if introduced_in_db(h, source_hashes) == db) for db in ("DB1", "DB2", "DB3")}
    with (args.outdir / "database_summary.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "database", "unique_sequences", "introduced_sequences", "total_aa", "min_aa", "mean_aa", "max_aa"])
        for db in ("DB1", "DB2", "DB3"):
            members = db_members[db]
            lengths = [len(seqs[h]) for h in members]
            w.writerow([sample_id, db, len(members), introduced_counts[db], sum(lengths), min(lengths) if lengths else 0, (sum(lengths)/len(lengths)) if lengths else 0, max(lengths) if lengths else 0])

    with (args.outdir / "source_summary.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "source", "raw_records", "unique_sequences"])
        for source, _ in sources:
            w.writerow([sample_id, source, raw_counts[source], len(source_hashes[source])])

    # Explicit orientation accounting for both raw RNA ORFs and non-redundant AA entities.
    raw_sense = sum(1 for recs in provenance.values() for r in recs if r.source == "RNAseq" and r.orientation == "sense")
    raw_reverse = sum(1 for recs in provenance.values() for r in recs if r.source == "RNAseq" and r.orientation == "reverse")
    unique_orient = defaultdict(int)
    for h in db1:
        unique_orient[rnaseq_orientation_label(provenance[h])] += 1
    with (args.outdir / "orientation_summary.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "metric", "orientation", "count"])
        w.writerow([sample_id, "raw_rnaseq_orfs", "sense", raw_sense])
        w.writerow([sample_id, "raw_rnaseq_orfs", "reverse", raw_reverse])
        for orientation in ("sense", "reverse", "sense+reverse"):
            w.writerow([sample_id, "unique_aa_sequences", orientation, unique_orient.get(orientation, 0)])

    with (args.outdir / "database_manifest.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "proteogenscope_version", "profile", "entity_type", "dataset_id", "patient_id", "sample_id",
            "tnm_label", "tnm_basis", "orf_mode", "orf_min_nt", "orf_reverse", "database", "definition",
            "fasta", "fasta_sha256", "unique_sequences", "rnaseq_orf_fasta", "rnaseq_orf_sha256",
            "canonical_fasta", "canonical_sha256", "isoforms_fasta", "isoforms_sha256",
        ])
        rnaseq_sha, canonical_sha, isoforms_sha = file_hash(args.rnaseq), file_hash(args.canonical), file_hash(args.isoforms)
        for db in ("DB1", "DB2", "DB3"):
            fp = fasta_paths[db]
            w.writerow([
                PROTEOGENSCOPE_VERSION, "proteoforms", ENTITY_TYPE, meta["dataset_id"], meta["patient_id"], sample_id,
                meta["tnm_label"], meta["tnm_basis"], args.orf_mode, args.orf_min_nt, args.orf_reverse.upper(), db,
                DB_DEFINITIONS[db], fp.name, file_hash(fp), len(db_members[db]), args.rnaseq.name, rnaseq_sha,
                args.canonical.name, canonical_sha, args.isoforms.name, isoforms_sha,
            ])

    print(f"ProteogenScope {PROTEOGENSCOPE_VERSION}")
    print(f"Profile: proteoforms | Sample: {sample_id} | Patient: {meta['patient_id']}")
    print(f"ORFs: mode={args.orf_mode} min_nt={args.orf_min_nt} reverse={args.orf_reverse.upper()}")
    print(f"RNA ORF orientation: sense={raw_sense} reverse={raw_reverse}")
    for db in ("DB1", "DB2", "DB3"):
        print(f"{db}: {len(db_members[db])} unique AA sequences ({introduced_counts[db]} introduced) -> {fasta_paths[db].name}")
    print(f"Outputs: {args.outdir.resolve()}")


if __name__ == "__main__":
    main()
