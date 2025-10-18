#!/usr/bin/env python3
"""Convert HYMET classified_sequences.tsv into a CAMI taxonomic profile."""

from __future__ import annotations

import csv
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
RANK_ALIAS = {
    "domain": "superkingdom",
    "kingdom": "superkingdom",
    "sk": "superkingdom",
    "k": "superkingdom",
    "phylum": "phylum",
    "p": "phylum",
    "class": "class",
    "c": "class",
    "order": "order",
    "o": "order",
    "family": "family",
    "f": "family",
    "genus": "genus",
    "g": "genus",
    "species": "species",
    "s": "species",
}


def _run(cmd: List[str], input_text: str | None = None) -> str:
    proc = subprocess.run(
        cmd,
        input=input_text,
        text=True,
        capture_output=True,
        check=True,
    )
    return proc.stdout


def parse_lineage(lineage: str) -> Dict[str, str]:
    out = {rank: "" for rank in RANKS}
    if not lineage:
        return out
    for part in lineage.split(";"):
        part = part.strip()
        if not part or ":" not in part:
            continue
        rk, name = part.split(":", 1)
        rk = RANK_ALIAS.get(rk.strip().lower(), rk.strip().lower())
        if rk in out:
            out[rk] = name.strip()
    return out


def batch_name2taxid(names: Iterable[str], taxdb: str) -> Dict[str, str]:
    names = sorted({n for n in names if n})
    if not names:
        return {}
    cmd = ["taxonkit", "name2taxid", "--data-dir", taxdb, "--show-rank"]
    input_text = "\n".join(names) + "\n"
    output = _run(cmd, input_text)
    mapping: Dict[str, str] = {}
    for line in output.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            mapping[parts[0]] = parts[1]
    return mapping


def batch_taxpath(taxids: Iterable[str], taxdb: str) -> Dict[str, Tuple[str, str]]:
    taxids = sorted({t for t in taxids if t})
    if not taxids:
        return {}
    cmd = [
        "taxonkit",
        "reformat",
        "--data-dir",
        taxdb,
        "-I",
        "1",
        "-f",
        "{d}|{p}|{c}|{o}|{f}|{g}|{s}",
        "-t",
    ]
    input_text = "\n".join(taxids) + "\n"
    output = _run(cmd, input_text)
    mapping: Dict[str, Tuple[str, str]] = {}
    for line in output.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 3:
            mapping[parts[0]] = (parts[1], parts[2])
    return mapping


def load_records(path: Path) -> List[Dict[str, str]]:
    records: List[Dict[str, str]] = []
    with path.open(encoding="utf-8", errors="ignore") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            parsed = parse_lineage(row.get("Lineage", ""))
            if any(parsed.values()):
                records.append(parsed)
    return records


def accumulate(records: List[Dict[str, str]], name2tid: Dict[str, str]) -> Tuple[Dict[str, Dict[str, int]], Dict[str, int], Iterable[str]]:
    counts = {rank: defaultdict(int) for rank in RANKS}
    totals = {rank: 0 for rank in RANKS}
    taxids_needed: set[str] = set()
    for parsed in records:
        for rank in RANKS:
            name = parsed.get(rank)
            if not name:
                continue
            tid = name2tid.get(name)
            if not tid:
                continue
            counts[rank][tid] += 1
            totals[rank] += 1
            taxids_needed.add(tid)
    return counts, totals, taxids_needed


def emit_cami(counts: Dict[str, Dict[str, int]], totals: Dict[str, int], taxid2path: Dict[str, Tuple[str, str]]) -> None:
    print("#CAMI Submission for Taxonomic Profiling")
    print("@Version:0.9.1 @Ranks:superkingdom|phylum|class|order|family|genus|species @SampleID:sample_0")
    print("@@TAXID RANK TAXPATH TAXPATHSN PERCENTAGE")
    for rank in RANKS:
        total = totals.get(rank, 0)
        if total <= 0:
            continue
        rank_counts = counts.get(rank, {})
        for tid, count in sorted(rank_counts.items(), key=lambda kv: kv[1], reverse=True):
            path = taxid2path.get(tid)
            if not path:
                continue
            names, ids = path
            pct = 100.0 * count / total
            print(f"{tid}\t{rank}\t{ids}\t{names}\t{pct:.6f}")


def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: hymet2cami.py <classified_sequences.tsv>", file=sys.stderr)
        sys.exit(1)

    input_path = Path(sys.argv[1]).resolve()
    if not input_path.is_file():
        print(f"Missing classified_sequences TSV: {input_path}", file=sys.stderr)
        sys.exit(1)

    taxdb = os.environ.get("TAXONKIT_DB", str(Path(__file__).resolve().parents[1] / "taxonomy_files"))
    print(f"[hymet2cami] using taxonomy DB at {taxdb}", file=sys.stderr)

    records = load_records(input_path)
    print(f"[hymet2cami] parsed {len(records)} lineages", file=sys.stderr)

    all_names = set()
    for parsed in records:
        for name in parsed.values():
            if name:
                all_names.add(name)

    print(f"[hymet2cami] converting {len(all_names)} unique taxon names", file=sys.stderr)
    name2tid = batch_name2taxid(all_names, taxdb)
    print(f"[hymet2cami] mapped {len(name2tid)} names to taxids", file=sys.stderr)

    counts, totals, needed_taxids = accumulate(records, name2tid)
    print(f"[hymet2cami] converting {len(needed_taxids)} taxids to paths", file=sys.stderr)
    taxid2path = batch_taxpath(needed_taxids, taxdb)

    emit_cami(counts, totals, taxid2path)
    print("[hymet2cami] done", file=sys.stderr)



if __name__ == "__main__":
    main()
