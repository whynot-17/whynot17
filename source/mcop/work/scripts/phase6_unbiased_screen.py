from __future__ import annotations

import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "phase6_unbiased_screen"
CACHE = WORK / "signature_cache"
WORK.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)
META_URL = "https://maayanlab.cloud/sigcom-lincs/metadata-api/signatures/find"

sys.path.insert(0, str(Path(__file__).resolve().parent))
import phase5_perturbation_reversal as p5  # noqa: E402
import phase3_module_decomposition as p3  # noqa: E402

# A priori approved, non-oncology screening universe. The list is intentionally
# broad across neurology, cardiovascular, metabolic, infectious, GI and other
# non-cancer indications; oncology drugs are retained only as internal controls.
APPROVED_NON_ONCOLOGY = sorted({
    "acetazolamide", "acyclovir", "albendazole", "allopurinol", "amantadine", "amitriptyline",
    "amiodarone", "amlodipine", "amoxicillin", "atenolol", "atorvastatin", "azithromycin",
    "baclofen", "budesonide", "buspirone", "caffeine", "captopril", "carbamazepine",
    "celecoxib", "chloroquine", "chlorpromazine", "cimetidine", "ciprofloxacin", "citalopram",
    "clomipramine", "clonazepam", "clopidogrel", "colchicine", "cyclosporine", "dapsone",
    "dexamethasone", "diazepam", "diclofenac", "digoxin", "diphenhydramine", "disulfiram",
    "doxycycline", "duloxetine", "entacapone", "erythromycin", "escitalopram", "estradiol",
    "fluoxetine", "fluvoxamine", "furosemide", "gabapentin", "gemfibrozil", "glibenclamide",
    "glimepiride", "haloperidol", "hydroxychloroquine", "hydroxyurea", "ibuprofen", "indomethacin",
    "itraconazole", "ivermectin", "ketoconazole", "lansoprazole", "leflunomide", "levodopa",
    "lithium", "losartan", "melatonin", "metformin", "methotrexate", "metoprolol", "miconazole",
    "minocycline", "montelukast", "morphine", "niclosamide", "nifedipine", "nitazoxanide",
    "omeprazole", "paroxetine", "phenformin", "pioglitazone", "pravastatin", "propranolol",
    "quetiapine", "rapamycin", "ribavirin", "rifampicin", "riluzole", "rosuvastatin", "salbutamol",
    "sertraline", "sildenafil", "simvastatin", "sitagliptin", "sodium phenylbutyrate", "spironolactone",
    "sulfasalazine", "tacrolimus", "teriflunomide", "thioridazine", "topiramate", "trazodone",
    "valproic acid", "valproate", "verapamil", "warfarin", "zidovudine", "zinc",
})

CONTROL_DRUGS = {"leflunomide", "teriflunomide", "bortezomib", "meldonium"}
CONTEXTS = ["HT29", "HCC515", "A549", "MCF7", "PC3"]
DOWNLOAD_MISSING = False
PRIORITY_NAMES = {
    "metformin", "chloroquine", "hydroxychloroquine", "niclosamide", "disulfiram", "valproic acid",
    "sodium phenylbutyrate", "dexamethasone", "rapamycin", "itraconazole", "propranolol", "melatonin",
    "teriflunomide", "leflunomide", "bortezomib", "meldonium",
}


def query_metadata(where: dict, limit: int = 1000) -> list[dict]:
    payload = {"filter": {"where": where, "limit": limit}}
    response = requests.post(META_URL, json=payload, timeout=25)
    response.raise_for_status()
    return response.json()


def get_ht29_metadata() -> list[dict]:
    path = WORK / "named_drug_metadata.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))

    # The unrestricted HT29 page is dominated by anonymous BRD IDs. Querying
    # by the curated drug name is slower but preserves interpretable drug
    # identities and avoids pretending an anonymous research compound is an
    # approved medicine.
    names = sorted(PRIORITY_NAMES)

    def one(name: str) -> list[dict]:
        try:
            records = query_metadata({"meta.pert_name": name}, 1000)
            return [r for r in records if r.get("meta", {}).get("cell_line") in CONTEXTS]
        except Exception:
            return []

    records: list[dict] = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(one, name) for name in names]
        for future in as_completed(futures):
            records.extend(future.result())
    # Keep one copy of each deposited signature.
    unique = {r.get("id"): r for r in records}
    records = list(unique.values())
    path.write_text(json.dumps(records), encoding="utf-8")
    return records


def choose_candidate_records(records: list[dict]) -> pd.DataFrame:
    rows = []
    for r in records:
        m = r.get("meta", {})
        drug = str(m.get("pert_name", "")).strip().lower()
        if not drug or m.get("pert_type") != "Chemical":
            continue
        # Candidate universe is the curated approved non-oncology list plus
        # three prespecified controls. We retain all HT29 times/doses.
        if drug not in APPROVED_NON_ONCOLOGY and drug not in CONTROL_DRUGS:
            continue
        rows.append({"drug": drug, "signature_id": m.get("local_id", r.get("id", "")), "cell_line": m.get("cell_line"), "pert_time": m.get("pert_time"), "pert_dose": m.get("pert_dose"), "persistent_id": m.get("persistent_id"), "tissue": m.get("tissue"), "disease": m.get("disease")})
    frame = pd.DataFrame(rows).drop_duplicates("signature_id")
    if frame.empty:
        return frame
    # Balanced, tractable perturbation panel: at most four signatures per
    # drug/context, while favoring 24 h and retaining dose/time diversity.
    frame["_time_priority"] = (frame["pert_time"].astype(str) != "24 h").astype(int)
    frame = frame.sort_values(["drug", "cell_line", "_time_priority", "pert_time", "pert_dose", "signature_id"])
    frame = frame.groupby(["drug", "cell_line"], group_keys=False).head(4).drop(columns=["_time_priority"])
    return frame.reset_index(drop=True)


def read_signature(meta: dict) -> pd.Series:
    url = meta["persistent_id"]
    filename = CACHE / Path(url).name
    if not filename.exists():
        if not DOWNLOAD_MISSING:
            raise FileNotFoundError(filename)
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        filename.write_bytes(response.content)
    frame = pd.read_csv(filename, sep="\t")
    frame.columns = [str(x).strip() for x in frame.columns]
    symbol_col, value_col = frame.columns[:2]
    frame[symbol_col] = frame[symbol_col].astype(str).str.upper().str.strip()
    frame[value_col] = pd.to_numeric(frame[value_col], errors="coerce")
    frame = frame.dropna(subset=[symbol_col, value_col])
    return frame.groupby(symbol_col)[value_col].median()


def load_disease_matrix() -> pd.DataFrame:
    path = ROOT / "phase5_perturbation_reversal" / "gene_delta_matrix_primary.csv"
    return pd.read_csv(path, index_col=0)


def load_state() -> pd.DataFrame:
    path = ROOT / "phase4_vulnerability_mapping" / "model_vulnerability_cooccurrence.csv"
    state = pd.read_csv(path).set_index("contrast_id")
    return state.reindex(p5.CONTRAST_IDS)


def top_genes(disease: pd.Series, n_each: int = 250) -> list[str]:
    d = disease.dropna()
    return list(dict.fromkeys(list(d.nlargest(min(n_each, len(d))).index) + list(d.nsmallest(min(n_each, len(d))).index)))


def score(disease: pd.Series, drug: pd.Series, genes: list[str]) -> tuple[float, float, int]:
    joined = pd.concat([disease.reindex(genes).rename("d"), drug.reindex(genes).rename("p")], axis=1).dropna()
    if len(joined) < 20:
        return np.nan, np.nan, len(joined)
    rho = float(spearmanr(joined["d"], joined["p"]).statistic)
    return (1 - rho) / 2, rho, len(joined)


def signature_worker(row: dict) -> tuple[str, pd.Series | None]:
    try:
        return row["signature_id"], read_signature(row)
    except Exception:
        return row["signature_id"], None


def build_subtype_diseases(matrix: pd.DataFrame, state: pd.DataFrame) -> dict[str, pd.Series]:
    diseases = {"global": matrix[p5.CONTRAST_IDS].median(axis=1)}
    rules = {
        "DHODH_subtype": state["salvage_low_DHODH_high"].fillna(False).astype(bool),
        "ERAD_subtype": state["UPR_low_ERAD_high"].fillna(False).astype(bool),
        "EMT_subtype": state["mesenchymal_high"].fillna(False).astype(bool),
    }
    for name, mask in rules.items():
        ids = [x for x in p5.CONTRAST_IDS if bool(mask.get(x, False))]
        diseases[name] = matrix[ids].median(axis=1)
    return diseases


def main() -> None:
    records = get_ht29_metadata()
    candidates = choose_candidate_records(records)
    candidates.to_csv(WORK / "candidate_signature_metadata_ht29.csv", index=False)
    if candidates.empty:
        raise SystemExit("No candidates found in the first HT29 metadata page")
    sigs: dict[str, pd.Series] = {}
    rows = candidates.to_dict("records")
    with ThreadPoolExecutor(max_workers=32) as pool:
        futures = [pool.submit(signature_worker, row) for row in rows]
        for future in as_completed(futures):
            sid, sig = future.result()
            if sig is not None:
                sigs[sid] = sig
    matrix = load_disease_matrix()
    state = load_state()
    diseases = build_subtype_diseases(matrix, state)
    global_genes = {name: top_genes(d) for name, d in diseases.items()}
    records_out = []
    for row in rows:
        sid = row["signature_id"]
        sig = sigs.get(sid)
        if sig is None:
            continue
        out = dict(row)
        for subtype, disease in diseases.items():
            s, rho, n = score(disease, sig, global_genes[subtype])
            out[f"rrs_{subtype}"] = s
            out[f"rho_{subtype}"] = rho
            out[f"n_{subtype}"] = n
        records_out.append(out)
    scored = pd.DataFrame(records_out)
    scored.to_csv(WORK / "per_signature_subtype_rrs.csv", index=False)
    grouped = scored.groupby("drug").agg(
        n_signatures=("signature_id", "count"),
        multi_context_RRS=("rrs_global", "median"),
        DHODH_subtype_RRS=("rrs_DHODH_subtype", "median"),
        ERAD_subtype_RRS=("rrs_ERAD_subtype", "median"),
        EMT_subtype_RRS=("rrs_EMT_subtype", "median"),
        multi_context_IQR=("rrs_global", lambda x: x.quantile(.75) - x.quantile(.25)),
        DHODH_IQR=("rrs_DHODH_subtype", lambda x: x.quantile(.75) - x.quantile(.25)),
        ERAD_IQR=("rrs_ERAD_subtype", lambda x: x.quantile(.75) - x.quantile(.25)),
    ).reset_index()
    ht29 = scored[scored["cell_line"] == "HT29"].groupby("drug")["rrs_global"].median().rename("HT29_global_RRS")
    context_stats = scored.groupby(["drug", "cell_line"])["rrs_global"].median().reset_index()
    context_summary = context_stats.groupby("drug")["rrs_global"].agg(
        n_contexts="count", context_RRS_SD="std", context_RRS_min="min", context_RRS_max="max"
    )
    grouped = grouped.merge(ht29, on="drug", how="left").merge(context_summary, on="drug", how="left")
    for subtype in ["DHODH_subtype", "ERAD_subtype", "EMT_subtype"]:
        grouped[f"{subtype}_selectivity"] = grouped[f"{subtype}_RRS"] - grouped["multi_context_RRS"]
    grouped["role"] = grouped["drug"].map({"leflunomide": "DHODH comparator", "teriflunomide": "DHODH comparator", "bortezomib": "proteostasis positive control", "meldonium": "failed appendix candidate"}).fillna("approved non-oncology panel")
    grouped["lead_gate"] = (
        (grouped["n_contexts"] >= 4)
        & (grouped["n_signatures"] >= 10)
        & (grouped["ERAD_subtype_RRS"] >= 0.60)
        & (grouped["ERAD_subtype_selectivity"] >= 0.10)
    )
    grouped["cross_context_confidence"] = np.where(
        (grouped["n_contexts"] >= 4) & (grouped["context_RRS_SD"] <= 0.05), "supported",
        np.where(grouped["n_contexts"] >= 3, "exploratory", "insufficient_contexts")
    )
    grouped.sort_values("ERAD_subtype_RRS", ascending=False).to_csv(WORK / "drug_subtype_rrs_ranked.csv", index=False)
    manifest = {"ht29_metadata_records": len(records), "candidate_metadata_records": len(candidates), "scored_signatures": len(scored), "unique_scored_drugs": int(scored["drug"].nunique()), "contexts": CONTEXTS, "note": "This is the first HT29 page of the LINCS metadata universe; approved status is a curated filter and requires formal FDA/ChEMBL verification before publication."}
    (WORK / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    print(grouped.sort_values("ERAD_subtype_RRS", ascending=False).head(20).to_string(index=False))


if __name__ == "__main__":
    main()
