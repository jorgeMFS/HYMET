#!/usr/bin/env python3
"""Convert MetaPhlAn4 profile outputs into CAMI format."""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Tuple

if __package__ is None or __package__ == "":  # pragma: no cover - CLI fallback
    sys.path.append(os.path.dirname(__file__))
    from common import RANKS, taxonkit_name2taxid, taxonkit_taxpath, write_cami_profile  # type: ignore
else:  # pragma: no cover
    from .common import RANKS, taxonkit_name2taxid, taxonkit_taxpath, write_cami_profile


def read_metaphlan(path: str) -> List[Tuple[str, float]]:
    rows: List[Tuple[str, float]] = []
    with open(path, "r") as handle:
        for raw in handle:
            raw = raw.strip()
            if not raw or raw.startswith("#"):
                continue
            parts = raw.split("\t")
            if len(parts) < 2:
                continue
            lineage = parts[0].strip()
            try:
                abundance = float(parts[1])
            except ValueError:
                continue
            rows.append((lineage, abundance))
    return rows


def lineage_to_ranked_names(lineage: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    components = lineage.split("|")
    for comp in components:
        if "__" not in comp:
            continue
        prefix, name = comp.split("__", 1)
        name = name.replace("_", " ").strip()
        prefix = prefix.lower()
        if prefix == "k":
            out["superkingdom"] = name
        elif prefix == "p":
            out["phylum"] = name
        elif prefix == "c":
            out["class"] = name
        elif prefix == "o":
            out["order"] = name
        elif prefix == "f":
            out["family"] = name
        elif prefix == "g":
            out["genus"] = name
        elif prefix == "s":
            out["species"] = name
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert MetaPhlAn4 profiles to CAMI format.")
    ap.add_argument("--input", required=True, help="MetaPhlAn profile TSV.")
    ap.add_argument("--out", required=True, help="Output CAMI TSV.")
    ap.add_argument("--sample-id", required=True, help="Sample identifier.")
    ap.add_argument("--tool", default="metaphlan4", help="Tool identifier.")
    ap.add_argument("--taxdb", default=os.environ.get("TAXONKIT_DB", ""), help="TaxonKit database directory.")
    args = ap.parse_args()

    records = read_metaphlan(args.input)
    if not records:
        write_cami_profile([], args.out, args.sample_id, args.tool, normalise=False)
        return

    final_names = []
    for lineage, _ in records:
        ranked = lineage_to_ranked_names(lineage)
        for rank in reversed(RANKS):
            if rank in ranked:
                final_names.append(ranked[rank])
                break

    name_to_taxid = taxonkit_name2taxid(final_names, args.taxdb)
    taxids = [taxid for taxid, _ in name_to_taxid.values()]
    taxid_to_paths = taxonkit_taxpath(taxids, args.taxdb)

    cami_rows = []
    for lineage, abundance in records:
        ranked = lineage_to_ranked_names(lineage)
        taxid = "NA"
        for rank in reversed(RANKS):
            nm = ranked.get(rank)
            if nm and nm in name_to_taxid:
                candidate, c_rank = name_to_taxid[nm]
                if candidate in taxid_to_paths:
                    taxid = candidate
                    ids_str, names_str = taxid_to_paths[candidate]
                    taxpath = ids_str.split("|")
                    names = names_str.split("|")
                    break
        else:
            taxpath = ["NA"] * len(RANKS)
            names = ["NA"] * len(RANKS)
        target_rank = None
        for rank in reversed(RANKS):
            if rank in ranked:
                target_rank = rank
                break
        cami_rows.append(
            {
                "taxid": taxid,
                "rank": target_rank or "species",
                "taxpath": taxpath,
                "taxpathsn": names,
                "percentage": abundance,
            }
        )

    write_cami_profile(cami_rows, args.out, args.sample_id, args.tool, normalise=True)


if __name__ == "__main__":
    main()
