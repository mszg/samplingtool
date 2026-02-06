"""Lightweight taxon name/ID validation against the NCBI taxonomy JSONL dump.

This module keeps a cached, lower‑cased index of scientific names -> taxids and a
simple taxid -> record map. It deliberately avoids building the full sampling
"tree" so that we can validate user input quickly before kicking off expensive
jobs.
"""

from __future__ import annotations

import json
from collections import defaultdict
from functools import lru_cache
from typing import Iterable, List, Dict, Any
import os


def _normalise_name(name: str) -> str:
    return name.strip().lower()


@lru_cache(maxsize=4)
def _load_index(taxonomy_path: str):
    if not taxonomy_path:
        raise FileNotFoundError("No taxonomy path provided")
    if not os.path.exists(taxonomy_path):
        raise FileNotFoundError(f"Taxonomy file not found: {taxonomy_path}")

    name_to_taxids: Dict[str, List[int]] = defaultdict(list)
    taxid_to_record: Dict[int, Dict[str, Any]] = {}

    with open(taxonomy_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            obj = json.loads(line)
            taxonomy = obj.get("taxonomy", {})
            tax_id = taxonomy.get("tax_id")
            sci_name = taxonomy.get("current_scientific_name", {}).get("name")
            rank = taxonomy.get("rank")

            if not tax_id or not sci_name:
                continue

            # canonical map for returning extra metadata
            taxid_to_record[tax_id] = {
                "tax_id": tax_id,
                "name": sci_name,
                "rank": rank,
            }

            name_to_taxids[_normalise_name(sci_name)].append(tax_id)

            # include known synonyms if present in the JSONL
            for syn in taxonomy.get("synonyms", []) or []:
                if not syn:
                    continue
                name_to_taxids[_normalise_name(syn)].append(tax_id)

    return name_to_taxids, taxid_to_record


def resolve_to_taxids(query: str | int, taxonomy_path: str) -> List[Dict[str, Any]]:
    """Return one or more taxonomy records that match the query.

    - Numeric input is treated as a taxid and validated for existence.
    - String input is matched case‑insensitively against scientific names and
      synonyms.
    - Returns a list of record dicts (tax_id, name, rank). Empty list means no
      match.
    """
    if query is None:
        return []

    name_to_taxids, taxid_to_record = _load_index(taxonomy_path)

    if isinstance(query, int) or (isinstance(query, str) and query.isdigit()):
        tid = int(query)
        rec = taxid_to_record.get(tid)
        return [rec] if rec else []

    key = _normalise_name(str(query))
    if not key:
        return []

    taxids = name_to_taxids.get(key, [])
    return [taxid_to_record[tid] for tid in taxids if tid in taxid_to_record]
