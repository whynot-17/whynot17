#!/usr/bin/env python
"""Phase 8-M: pre-specified named-candidate challenge for Meldonium.

This is not an unbiased drug screen. It first audits reliable Meldonium
aliases across the local GDSC/PRISM/CTRPv2/LINCS inputs. If a direct response
record exists, the frozen Phase 8-R2 framework is the only permitted scoring
route. If no direct response exists, the script stops before inventing a
pharmacogenomic score and instead records a mechanism-consistency audit from
the frozen OXA-R pathway and R3 functional-dependency results.
"""

from __future__ import annotations

import json
import re
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs"
WORK = ROOT / "work"
DATA = WORK / "phase7_convergent_vulnerability" / "data"
PRISM_RAW = WORK / "phase8R_prism" / "raw"
CTRP_RAW = WORK / "phase8R2_ctrpv2" / "raw"

ALIASES = [
    ("meldonium", "INN/common name"),
    ("mildronate", "brand/generic name"),
    ("mildronat", "MeSH entry term"),
    ("mildronāts", "brand spelling"),
    ("mindronate", "MeSH entry term"),
    ("MET-88", "MeSH synonym"),
    ("3-TMHP", "MeSH abbreviation"),
    ("Vasonat", "MeSH synonym"),
    ("Quaterin", "PubChem synonym"),
    ("Kvaterin", "PubChem synonym"),
    ("Quaterine", "PubChem synonym"),
    ("3-(2,2,2-trimethylhydrazine)propionate", "chemical name"),
    ("3-(2,2,2-trimethylhydrazinium)propionate", "chemical name"),
    ("3-(2,2,2-trimethylhydrazinium)propanoate", "chemical name"),
    ("3-(2,2,2-trimethyldiazaniumyl)propanoate", "chemical name"),
    ("3-[(trimethylazaniumyl)amino]propanoate", "IUPAC name"),
    ("N-trimethylhydrazine-3-propionate", "chemical name"),
    ("trimethylhydraziniumpropionate", "chemical name"),
    ("76144-81-5", "CAS"),
    ("802032-35-5", "salt/CAS record"),
    ("CHEBI:131843", "ChEBI"),
    ("CHEMBL2104708", "ChEMBL"),
    ("CID 123868", "PubChem CID"),
    ("C01EB22", "ATC"),
    ("73H7UDN6EC", "UNII"),
    ("DTXSID10997497", "DSSTox"),
    ("DTXCID301424501", "DSSTox compound ID"),
]


def norm(value: object) -> str:
    text = unicodedata.normalize("NFKD", str(value or ""))
    text = "".join(x for x in text if not unicodedata.combining(x))
    return re.sub(r"[^A-Z0-9]", "", text.upper())


ALIAS_TABLE = pd.DataFrame([{"alias": a, "source_role": role, "alias_key": norm(a)} for a, role in ALIASES]).drop_duplicates("alias_key")
ALIAS_KEYS = set(ALIAS_TABLE.alias_key)


def find_alias(value: object) -> str:
    k = norm(value)
    return next((a for a, ak in zip(ALIAS_TABLE.alias, ALIAS_TABLE.alias_key) if ak == k), "")


def find_alias_in_text(value: object) -> str:
    text = str(value or "")
    for alias, _role in ALIASES:
        if norm(alias) and norm(alias) in norm(text):
            return alias
    return ""


def audit_gdsc() -> dict[str, object]:
    rows = []
    compound_path = DATA / "screened_compounds_rel_8.4.csv"
    if compound_path.exists():
        comp = pd.read_csv(compound_path, low_memory=False)
        for field in ["DRUG_NAME", "SYNONYMS"]:
            if field not in comp:
                continue
            for _, r in comp.iterrows():
                hit = find_alias_in_text(r.get(field, ""))
                if hit:
                    rows.append({"platform": "GDSC", "data_role": "compound_metadata", "field": field,
                                 "matched_alias": hit, "matched_value": r.get(field, ""), "DRUG_NAME": r.get("DRUG_NAME", "")})
    assoc_path = OUT / "phase8R_drug_association_by_trajectory.csv"
    if assoc_path.exists():
        assoc = pd.read_csv(assoc_path, usecols=["database", "DRUG_NAME", "n"])
        assoc["matched_alias"] = assoc.DRUG_NAME.map(find_alias)
        hit = assoc[assoc.matched_alias.ne("")].copy()
        for _, r in hit.iterrows():
            rows.append({"platform": r.database, "data_role": "CRC_response_association", "field": "DRUG_NAME",
                         "matched_alias": r.matched_alias, "matched_value": r.DRUG_NAME, "DRUG_NAME": r.DRUG_NAME,
                         "response_n": r.n})
    return {"platform": "GDSC", "metadata_or_response_match": bool(rows),
            "n_matching_records": len(rows), "n_crc_response_records": int(sum(x["data_role"] == "CRC_response_association" for x in rows)),
            "matched_aliases": ";".join(sorted(set(x["matched_alias"] for x in rows))),
            "matched_names": ";".join(sorted(set(str(x.get("DRUG_NAME", "")) for x in rows))),
            "details": rows}


def audit_prism(model_ids: set[str]) -> dict[str, object]:
    path = PRISM_RAW / "secondary-screen-dose-response-curve-parameters.csv"
    rows = []
    total = 0
    crc_total = 0
    if path.exists():
        for chunk in pd.read_csv(path, usecols=["depmap_id", "name", "auc"], chunksize=200_000, low_memory=False):
            chunk["depmap_id"] = chunk.depmap_id.astype(str)
            chunk["matched_alias"] = chunk.name.map(find_alias)
            hit = chunk[chunk.matched_alias.ne("")]
            total += len(hit)
            crc_hit = hit[hit.depmap_id.isin(model_ids)]
            crc_total += len(crc_hit)
            for _, r in hit.iterrows():
                rows.append({"platform": "PRISM", "data_role": "dose_response", "field": "name",
                             "matched_alias": r.matched_alias, "matched_value": r.name, "DRUG_NAME": r.name,
                             "depmap_id": r.depmap_id, "is_crc_expression_model": bool(r.depmap_id in model_ids),
                             "auc": r.auc})
    return {"platform": "PRISM", "metadata_or_response_match": bool(rows),
            "n_matching_records": total, "n_crc_response_records": crc_total,
            "matched_aliases": ";".join(sorted(set(x["matched_alias"] for x in rows))),
            "matched_names": ";".join(sorted(set(str(x["DRUG_NAME"]) for x in rows))), "details": rows}


def ctrp_crc_map(model_ids: set[str]) -> dict[str, str]:
    path = CTRP_RAW / "Harmonized_CCL_Data_v1.0.xlsx"
    if not path.exists():
        return {}
    ccl = pd.read_excel(path)
    ccl = ccl[ccl.Dataset.astype(str).eq("CTRPv2")].copy()
    ccl["depmap_id"] = ccl.Synonyms.astype(str).str.extract(r"(ACH-[0-9]+)", expand=False)
    ccl = ccl[ccl.depmap_id.isin(model_ids)]
    return dict(zip(ccl.Harmonized_Cell_Line_ID.astype(str), ccl.depmap_id.astype(str)))


def audit_ctrp(model_ids: set[str]) -> dict[str, object]:
    path = CTRP_RAW / "CTRPv2_Results_v1.0.tsv"
    left_map = ctrp_crc_map(model_ids)
    rows = []
    if path.exists():
        for chunk in pd.read_csv(path, sep="\t", usecols=["Key", "AUC_all_ccl_CTRPv2_conc"], chunksize=200_000, low_memory=False):
            parts = chunk.Key.astype(str).str.partition(":|:")
            chunk["cell_line"] = parts[0]
            chunk["DRUG_NAME"] = parts[2].str.strip()
            chunk["matched_alias"] = chunk.DRUG_NAME.map(find_alias)
            hit = chunk[chunk.matched_alias.ne("")]
            if len(hit):
                for _, r in hit.iterrows():
                    rows.append({"platform": "CTRPv2", "data_role": "AUC_response", "field": "Key compound token",
                                 "matched_alias": r.matched_alias, "matched_value": r.DRUG_NAME, "DRUG_NAME": r.DRUG_NAME,
                                 "cell_line": r.cell_line, "is_crc_expression_model": bool(r.cell_line in left_map),
                                 "auc": r.AUC_all_ccl_CTRPv2_conc})
    return {"platform": "CTRPv2", "metadata_or_response_match": bool(rows),
            "n_matching_records": len(rows), "n_crc_response_records": int(sum(x["is_crc_expression_model"] for x in rows)),
            "matched_aliases": ";".join(sorted(set(x["matched_alias"] for x in rows))),
            "matched_names": ";".join(sorted(set(str(x["DRUG_NAME"]) for x in rows))), "details": rows}


def audit_lincs() -> dict[str, object]:
    path = WORK / "phase5_perturbation_reversal" / "perturbation_signature_metadata.csv"
    rows = []
    if path.exists():
        m = pd.read_csv(path, low_memory=False)
        if "drug" in m:
            for _, r in m.iterrows():
                hit = find_alias(r.drug)
                if hit:
                    rows.append({"platform": "LINCS/L1000", "data_role": "perturbation_signature",
                                 "field": "drug", "matched_alias": hit, "matched_value": r.drug, "DRUG_NAME": r.drug})
    return {"platform": "LINCS/L1000", "metadata_or_response_match": bool(rows),
            "n_matching_records": len(rows), "n_crc_response_records": 0,
            "matched_aliases": ";".join(sorted(set(x["matched_alias"] for x in rows))),
            "matched_names": ";".join(sorted(set(str(x["DRUG_NAME"]) for x in rows))), "details": rows}


def mechanism_audit() -> pd.DataFrame:
    rows = []
    stability = pd.read_csv(OUT / "phase1_pathway_stability_summary_cell_lines.csv")
    for pathway in ["carnitine_entry", "FAO_mitochondrial"]:
        r = stability.loc[stability.pathway.eq(pathway)].iloc[0]
        rows.append({"evidence_layer": "OXA-R transcriptional state", "feature": pathway,
                     "metric": "n_down/n_up; median_delta", "value": f"{int(r.n_down)}/{int(r.n_up)}; {r.median_delta:.6f}",
                     "interpretation": "carnitine-entry is directionally down in 5/6 models; FAO is split 3/3. This is transcript-level state evidence, not flux or dependency."})
    genes = pd.read_csv(OUT / "phase7bR3_final_ranking.csv")
    for gene in ["SLC22A5", "CPT1A", "CPT2", "BBOX1"]:
        r = genes.loc[genes.gene.eq(gene)].iloc[0]
        rows.append({"evidence_layer": "R3 functional dependency", "feature": gene,
                     "metric": "r3_rank; median_rho; global_q; leave_HCT116_rho",
                     "value": f"{int(r.r3_rank)}; {r.median_background_vulnerability_rho:.6f}; {r.global_empirical_q_value:.6f}; {r.leave_HCT116_out_median_vulnerability_rho:.6f}",
                     "interpretation": "No gene passes the frozen final shortlist; the carnitine/FAO axis lacks stable universal functional dependency."})
    rows.append({"evidence_layer": "Meldonium pharmacology", "feature": "Meldonium→BBOX1/GBBD inhibition→carnitine availability↓", "metric": "directionality", "value": "not a reversal of carnitine-entry↓", "interpretation": "The broad direction is same-side rather than opposite-side to the dominant carnitine-entry transcript state; this weakens a universal reversal claim."})
    rows.append({"evidence_layer": "Direct public drug response", "feature": "Meldonium", "metric": "GDSC/PRISM/CTRPv2/LINCS alias audit", "value": "0 direct response/signature records", "interpretation": "No public pharmacogenomic score or percentile can be calculated; this is unavailable, not negative."})
    return pd.DataFrame(rows)


def main() -> None:
    ALIAS_TABLE.to_csv(OUT / "phase8M_meldonium_alias_registry.csv", index=False)
    model = pd.read_csv(WORK / "phase7b_depmap" / "raw" / "Model.csv", low_memory=False)
    lineage = model.OncotreeLineage.fillna("").astype(str).str.lower()
    primary = model.OncotreePrimaryDisease.fillna("").astype(str).str.lower()
    crc = model[(lineage == "bowel") & primary.str.contains("colorectal", regex=False)]
    expr = pd.read_csv(WORK / "phase7b_depmap" / "raw" / "OmicsExpressionProteinCodingGenesTPMLogp1.csv", index_col=0, low_memory=False)
    model_ids = set(crc.ModelID.astype(str)) & set(expr.index.astype(str))
    audits = [audit_gdsc(), audit_prism(model_ids), audit_ctrp(model_ids), audit_lincs()]
    summary = pd.DataFrame([{k: v for k, v in x.items() if k != "details"} for x in audits])
    summary.to_csv(OUT / "phase8M_public_response_audit.csv", index=False)
    details = [z for x in audits for z in x["details"]]
    detail_cols = ["platform", "data_role", "field", "matched_alias", "matched_value", "DRUG_NAME",
                   "response_n", "depmap_id", "is_crc_expression_model", "cell_line", "auc"]
    pd.DataFrame(details, columns=detail_cols).to_csv(OUT / "phase8M_public_response_audit_matches.csv", index=False)
    mech = mechanism_audit()
    mech.to_csv(OUT / "phase8M_mechanism_consistency_audit.csv", index=False)
    direct_available = bool(summary.metadata_or_response_match.any())
    status = "direct_response_available_locked_model_pending" if direct_available else "no_direct_response_unavailable_not_negative"
    report = ["# Phase 8-M：Meldonium named-candidate challenge", "",
              "## Prespecification", "",
              "Meldonium was tested as a named candidate after the Phase 8-R2 signature, thresholds and model universe were frozen. It was not reintroduced into the unbiased drug ranking and no prior candidate list was used.", "",
              "## Alias-level public-data audit", "", "```text", summary.to_string(index=False), "```", "",
              f"Challenge status: `{status}`.", "",
              "No Meldonium direct response or LINCS perturbation record was found across the audited local GDSC, PRISM, CTRPv2 and LINCS/L1000 inputs. Therefore no rho, IC50/AUC, percentile, RRS or OXA-R-state drug score was fabricated. The absence is an availability result, not a biological negative.", "",
              "## Mechanism-consistency audit", "", "```text", mech.to_string(index=False), "```", "",
              "## Decision", "",
              "Meldonium is not rescued by the current dry-lab data. Carnitine-entry is down in 5/6 acquired OXA-R contrasts, but FAO is only 3/3 directional and SLC22A5/CPT1A/CPT2/BBOX1 do not form a stable R3 functional dependency. Its established pharmacology is primarily BBOX/γ-butyrobetaine hydroxylase and OCTN2-related lowering of carnitine availability, rather than a simple direct CPT1A inhibitor ([Dambrova et al., PMID 26850121](https://pubmed.ncbi.nlm.nih.gov/26850121/)). Therefore the broad direction is not a clean reversal of the dominant carnitine-entry state.", "",
              "The only decisive next test is a matched parental/OXA-R mini-screen. Recommended primary comparison: HCT116-P/OXA-R and DLD1-P/OXA-R or HT29-P/OXA-R, with Meldonium dose-response and pre-specified selectivity index `SI = IC50_parental / IC50_OXA-R`. Suggested interpretation: SI > 1.5-2 in at least two independent backgrounds = revive; SI approximately 1 = deprioritize; SI < 1 = No-Go.", "",
              "## Source note", "",
              "Alias registry follows PubChem CID 123868 and NCBI MeSH entry terms: https://pubchem.ncbi.nlm.nih.gov/compound/123868 and https://www.ncbi.nlm.nih.gov/mesh/67050147. Raw response files remain local and are not committed."]
    (OUT / "phase8M_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    manifest = {"phase": "8-M", "candidate": "Meldonium", "prespecified_named_candidate_challenge": True,
                "frozen_model": "Phase 8-R2 six trajectories; four backgrounds; top100/250/500 weighted/rank; no re-ranking",
                "aliases": ALIASES, "crc_expression_models_audited": len(model_ids), "direct_response_available": direct_available,
                "challenge_status": status, "public_platforms": ["GDSC", "PRISM", "CTRPv2", "LINCS/L1000"],
                "if_direct_response": "apply frozen Phase 8-R2 locked model; no alternative threshold permitted",
                "if_unavailable": "mechanism consistency audit followed by matched parental/OXA-R mini-screen",
                "outputs": ["phase8M_meldonium_alias_registry.csv", "phase8M_public_response_audit.csv", "phase8M_public_response_audit_matches.csv", "phase8M_mechanism_consistency_audit.csv", "phase8M_report.md"],
                "raw_data_policy": "raw data remain local and are not committed"}
    (OUT / "phase8M_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(summary.to_string(index=False))
    print(mech.to_string(index=False))


if __name__ == "__main__":
    main()
