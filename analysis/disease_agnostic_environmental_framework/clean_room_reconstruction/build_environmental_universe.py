"""Clean-room construction of the environmental chemical universe.

This module intentionally contains no disease, outcome, gene-set, or
association imports.  It reads only the CTD chemical vocabulary and the
frozen chemical-classification rules.
"""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import pandas as pd


def canonical_id(value: object) -> str:
    text = str(value).strip()
    return text[5:] if text.startswith("MESH:") else text


def read_ctd(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("# Fields:"):
                header = next(handle).lstrip("# ").rstrip("\r\n").split("\t")
                return pd.read_csv(handle, sep="\t", names=header, dtype=str,
                                    keep_default_na=False, comment="#", low_memory=False)
    raise ValueError(f"CTD field header not found: {path}")


def split_pipe(value: object) -> list[str]:
    text = str(value or "")
    return [item for item in text.split("|") if item]


def descendant(path: str, prefix: str) -> bool:
    return path == prefix or path.startswith(prefix + ".")


def classify(chemicals: pd.DataFrame, rules_path: Path,
             drugcentral_path: Path | None = None,
             pah_formula_path: Path | None = None) -> pd.DataFrame:
    rules = json.loads(rules_path.read_text(encoding="utf-8"))
    core_rules = rules["core_rules"]
    drug_regex = re.compile(rules.get("drug_exclusion_regex", r"(?!)"))

    drug_cas: set[str] = set()
    drug_inn: set[str] = set()
    if drugcentral_path and drugcentral_path.exists():
        drug = pd.read_csv(drugcentral_path, sep="\t", dtype=str, keep_default_na=False)
        if "CAS_RN" in drug.columns:
            drug_cas = {str(x).strip().lower() for x in drug["CAS_RN"] if str(x).strip()}
        if "INN" in drug.columns:
            drug_inn = {re.sub(r"[^a-z0-9]+", "", str(x).lower()) for x in drug["INN"] if str(x).strip()}

    formulas = {}
    if pah_formula_path and pah_formula_path.exists():
        formulas = json.loads(pah_formula_path.read_text(encoding="utf-8"))

    records: list[dict[str, object]] = []
    for row in chemicals.to_dict("records"):
        paths = []
        for column in ("TreeNumbers", "ParentTreeNumbers"):
            paths.extend(value.split("/", 1)[0] for value in split_pipe(row.get(column, "")))
        classes: list[str] = []
        prefixes: list[str] = []
        for category, category_prefixes in core_rules.items():
            hits = [prefix for prefix in category_prefixes if any(descendant(path, prefix) for path in paths)]
            if hits:
                classes.append(category)
                prefixes.extend(hits)

        drug_text = " ".join(str(row.get(column, "")) for column in
                              ("ChemicalName", "Definition", "MESHSynonyms", "CTDCuratedSynonyms"))
        semantic_drug = bool(drug_regex.search(drug_text))
        cas = str(row.get("CasRN", "")).strip().lower()
        normalized_name = re.sub(r"[^a-z0-9]+", "", str(row.get("ChemicalName", "")).lower())
        drugcentral_hit = bool((cas and cas in drug_cas) or
                               (normalized_name and normalized_name in drug_inn))

        cid = re.sub(r"^CID:", "", str(row.get("PubChemCID", "")).strip())
        formula = str(formulas.get(cid, ""))
        formula_valid = bool(re.fullmatch(r"C[0-9]+H[0-9]+", formula)) if pah_formula_path else pd.NA
        if semantic_drug or drugcentral_hit:
            classes, prefixes = [], []
        if "pahs" in classes and pah_formula_path and not formula_valid:
            classes = [x for x in classes if x != "pahs"]
            prefixes = [x for x in prefixes if x != "D02.455.426.559.847"]

        records.append({
            "chemical_id": canonical_id(row.get("ChemicalID", "")),
            "canonical_name": row.get("ChemicalName", ""),
            "CAS": row.get("CasRN", ""),
            "PubChemCID": row.get("PubChemCID", ""),
            "DTXSID": row.get("DTXSID", ""),
            "InChIKey": row.get("InChIKey", ""),
            "chemical_class": ";".join(sorted(set(classes))),
            "classification_rule_prefixes": ";".join(sorted(set(prefixes))),
            "eligible_as_environmental_exposure": bool(classes),
            "drug_like_exclusion": semantic_drug,
            "drugcentral_match": bool(drugcentral_hit),
            "pah_formula": formula,
            "pah_structure_valid": formula_valid,
            "classification_exclusion_reason": (
                "drugcentral_cas_or_inn" if drugcentral_hit else
                "drug_semantic_regex" if semantic_drug else ""
            ),
        })

    result = pd.DataFrame(records).drop_duplicates("chemical_id", keep="first")
    if result["chemical_id"].duplicated().any():
        raise AssertionError("Duplicate chemical IDs in clean-room universe")
    return result
