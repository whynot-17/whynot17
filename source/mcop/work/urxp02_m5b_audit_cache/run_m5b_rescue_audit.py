from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


REPO = Path(r"C:/Users/21634/Documents/Codex/2026-08-22/non/work/whynot17")
CACHE = REPO / "work/urxp02_m5b_audit_cache"
OUT = REPO / "analysis/disease_agnostic_environmental_framework/urxp02_m5b_single_cell_rescue_audit"

SOURCE_URLS = {
    "GSE182416": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182416",
    "GSE202109": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202109",
    "GSE183276": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183276",
    "GSE183277": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183277",
    "GSE207784": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207784",
    "GSE207784_PMC": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9613617/",
    "GSE189795": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189795",
    "GSE165824": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165824",
    "SCP1265": "https://singlecell.broadinstitute.org/single_cell/study/SCP1265/deep-learning-enables-genetic-analysis-of-the-human-thoracic-aorta",
    "KIDNEY_ATLAS": "https://www.nature.com/articles/s41586-023-05769-3",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_soft(path: Path) -> list[dict[str, str | list[str]]]:
    samples: list[dict[str, str | list[str]]] = []
    current: dict[str, str | list[str]] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    samples.append(current)
                current = {"accession": line.split("=", 1)[1].strip(), "characteristics": []}
            elif current is None:
                continue
            elif line.startswith("!Sample_title = "):
                current["title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_source_name_ch1 = "):
                current["source"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1 = "):
                current.setdefault("characteristics", []).append(line.split("=", 1)[1].strip())
    if current:
        samples.append(current)
    for sample in samples:
        chars = sample.get("characteristics", [])
        parsed: dict[str, str] = {}
        for item in chars if isinstance(chars, list) else []:
            if ":" in item:
                key, value = item.split(":", 1)
                parsed[key.strip().lower()] = value.strip()
        sample["parsed"] = parsed
    return samples


def soft_counts(samples: list[dict[str, str | list[str]]], sex_keys=("sex", "gender")) -> dict[str, object]:
    sex = []
    ages = []
    for s in samples:
        parsed = s.get("parsed", {})
        if not isinstance(parsed, dict):
            parsed = {}
        value = next((parsed[k].lower() for k in sex_keys if k in parsed), "")
        if value.startswith("male"):
            sex.append("M")
        elif value.startswith("female"):
            sex.append("F")
        else:
            sex.append("unknown")
        ages.append(any(k == "age" for k in parsed))
    return {
        "n": len(samples),
        "male": sum(v == "M" for v in sex),
        "female": sum(v == "F" for v in sex),
        "unknown": sum(v == "unknown" for v in sex),
        "age": sum(ages) == len(samples) if samples else False,
    }


def parse_kidney_metadata(path: Path) -> dict[str, object]:
    """Count donor/sample and healthy-reference subset from deposited metadata.

    The deposited tables have one unlabelled cell barcode field before the
    published header, so raw positions are used deliberately: patient/donor=24,
    condition.long=20, sex=29 (zero-based in the data row).  ``specimen`` at
    position 19 can split a donor into multiple specimen identifiers, whereas
    ``patient`` is the donor-level field used for the atlas donor counts.
    """
    donor_i, condition_i, sex_i = 24, 20, 29
    total_cells = 0
    all_donors: dict[str, str] = {}
    sex_cells: Counter[str] = Counter()
    ref_cells = 0
    ref_donors: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        for row in reader:
            if len(row) <= sex_i:
                continue
            total_cells += 1
            donor = row[donor_i].strip().strip('"')
            condition = row[condition_i].strip().strip('"')
            sex = row[sex_i].strip().strip('"').upper()
            if sex in {"M", "F"}:
                all_donors[donor] = sex
                sex_cells[sex] += 1
            if condition == "Normal Reference":
                ref_cells += 1
                ref_donors[donor] = sex
    return {
        "total_cells": total_cells,
        "total_donors": len(all_donors),
        "total_male": sum(v == "M" for v in all_donors.values()),
        "total_female": sum(v == "F" for v in all_donors.values()),
        "sex_cells": dict(sex_cells),
        "ref_cells": ref_cells,
        "ref_donors": len(ref_donors),
        "ref_male": sum(v == "M" for v in ref_donors.values()),
        "ref_female": sum(v == "F" for v in ref_donors.values()),
    }


def inspect_h5ad_obs(path: Path) -> dict[str, object]:
    """Read only h5ad dimensions and obs metadata; never access expression X."""
    import anndata as ad

    obj = ad.read_h5ad(path, backed="r")
    obs = obj.obs
    donor_col = "donor_id" if "donor_id" in obs.columns else None
    sex_col = "sex" if "sex" in obs.columns else None
    donor_sex: dict[str, str] = {}
    if donor_col and sex_col:
        for donor, sex in obs[[donor_col, sex_col]].drop_duplicates().itertuples(index=False, name=None):
            donor_sex[str(donor)] = str(sex).lower()
    result = {
        "n_cells": int(obj.n_obs),
        "n_genes": int(obj.n_vars),
        "donors": len(donor_sex),
        "male": sum(v == "male" for v in donor_sex.values()),
        "female": sum(v == "female" for v in donor_sex.values()),
        "donor_sex": donor_sex,
        "sex_column": sex_col,
        "donor_column": donor_col,
        "cell_type_annotations": "cell_type" in obs.columns or "cell_type__ontology_label" in obs.columns,
        "normal_label": "normal" in set(obs["disease__ontology_label"].astype(str)) if "disease__ontology_label" in obs.columns else False,
    }
    obj.file.close()
    return result


def criterion(male: int | str, female: int | str) -> tuple[str, str]:
    if not isinstance(male, int) or not isinstance(female, int):
        return "unknown", "unknown"
    preferred = male >= 3 and female >= 3
    return ("yes" if preferred else "no", "no" if preferred else "yes")


def row(**kwargs: object) -> dict[str, object]:
    return kwargs


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gse182416 = parse_soft(CACHE / "GSE182416_family.soft.gz")
    gse202109 = parse_soft(CACHE / "GSE202109_family.soft.gz")
    gse207784 = parse_soft(CACHE / "GSE207784_family.soft.gz")
    gse189795 = parse_soft(CACHE / "GSE189795_family.soft.gz")
    gse165824 = parse_soft(CACHE / "GSE165824_family.soft.gz")
    k276 = parse_kidney_metadata(CACHE / "GSE183276_metadata.txt.gz")
    k277 = parse_kidney_metadata(CACHE / "GSE183277_metadata.txt.gz")
    h165 = inspect_h5ad_obs(CACHE / "GSE165824_ascending_descending_human_aorta_v1.h5ad")
    c182 = soft_counts(gse182416)
    c202 = soft_counts(gse202109)
    c207 = soft_counts(gse207784)
    c189 = soft_counts(gse189795)
    c165 = soft_counts(gse165824)

    # Explicit control subsets are derived from the GEO sample-level disease field.
    controls207 = [s for s in gse207784 if str(s.get("parsed", {}).get("disease", "")).lower() == "control"]
    controls189 = [s for s in gse189795 if "control" in str(s.get("title", "")).lower()]
    c207_control = soft_counts(controls207)
    c189_control = soft_counts(controls189)

    rows = [
        row(dataset_id="GSE182416", dataset_alias="normal thyroid scRNA", tissue="thyroid", healthy_or_control="not healthy-only: normal-adjacent samples with PTC/follicular adenoma pathology", technology="10x scRNA-seq", n_donors_total=c182["n"], n_male=c182["male"], n_female=c182["female"], n_cells=54726, donor_id_available="sample/GSM-level only", sex_metadata_available="yes", raw_counts_available="yes: supplementary raw-count TXT", processed_matrix_available="yes", cell_type_annotations_available="yes", age_available="yes", eligible_for_pseudobulk_sex_DE="no: 0 male donors", eligible_for_cell_composition_analysis="no: 0 male donors", analysis_subset="all 7 samples", preferred_minimum_met="no", low_power_flag="yes", notes="GEO summary says 54,762 cells; supplementary filename says 54,726. All 7 samples are female; pathology is PTC or follicular adenoma, so this is not a healthy sex-balanced thyroid reference.", source_url=SOURCE_URLS["GSE182416"]),
        row(dataset_id="GSE202109", dataset_alias="healthy living-donor kidney", tissue="kidney", healthy_or_control="healthy living donor kidney", technology="10x scRNA-seq 3' v3", n_donors_total=c202["n"], n_male=c202["male"], n_female=c202["female"], n_cells=27677, donor_id_available="yes: donor-level sample metadata", sex_metadata_available="yes", raw_counts_available="yes: SRA raw reads and supplementary count files", processed_matrix_available="yes", cell_type_annotations_available="yes", age_available="yes", eligible_for_pseudobulk_sex_DE="yes", eligible_for_cell_composition_analysis="yes", analysis_subset="all 19 living donors", preferred_minimum_met="yes", low_power_flag="no", notes="9 male / 10 female donors; healthy cohort and preferred >=3 per sex criterion.", source_url=SOURCE_URLS["GSE202109"]),
        row(dataset_id="GSE183276", dataset_alias="Human Kidney Cell Atlas scCv3", tissue="kidney", healthy_or_control="mixed disease atlas with identifiable Normal Reference subset", technology="10x scRNA-seq", n_donors_total=k276["total_donors"], n_male=k276["total_male"], n_female=k276["total_female"], n_cells=k276["total_cells"], donor_id_available="yes: specimen/donor field", sex_metadata_available="yes", raw_counts_available="yes: deposited RDS count matrix", processed_matrix_available="yes", cell_type_annotations_available="yes", age_available="unknown/not in deposited metadata table", eligible_for_pseudobulk_sex_DE="yes for Normal Reference subset: 7M/13F", eligible_for_cell_composition_analysis="yes for Normal Reference subset", analysis_subset="Normal Reference: 20 donors (7M/13F), 21,650 cells", preferred_minimum_met="yes for Normal Reference", low_power_flag="no for healthy subset", notes="Total atlas: 45 donors (21M/24F), 109,741 cells. Disease samples are not to be pooled with the healthy reference for sex DE.", source_url=SOURCE_URLS["GSE183276"]),
        row(dataset_id="GSE183277", dataset_alias="Human Kidney Cell Atlas snCv3", tissue="kidney", healthy_or_control="mixed disease atlas with identifiable Normal Reference subset", technology="10x snRNA-seq", n_donors_total=k277["total_donors"], n_male=k277["total_male"], n_female=k277["total_female"], n_cells=k277["total_cells"], donor_id_available="yes: specimen/donor field", sex_metadata_available="yes", raw_counts_available="yes: deposited RDS count matrix", processed_matrix_available="yes", cell_type_annotations_available="yes", age_available="unknown/not in deposited metadata table", eligible_for_pseudobulk_sex_DE="yes for Normal Reference subset: 8M/10F", eligible_for_cell_composition_analysis="yes for Normal Reference subset", analysis_subset="Normal Reference: 18 donors (8M/10F), 88,460 cells", preferred_minimum_met="yes for Normal Reference", low_power_flag="no for healthy subset", notes="Total atlas: 36 donors (19M/17F), 203,702 cells. Disease samples are not to be pooled with the healthy reference for sex DE.", source_url=SOURCE_URLS["GSE183277"]),
        row(dataset_id="GSE207784", dataset_alias="aneurysmal human ascending aorta snRNA", tissue="ascending aorta", healthy_or_control="mixed: 7 controls and 6 aneurysm samples", technology="snRNA-seq", n_donors_total=c207["n"], n_male=c207["male"], n_female=c207["female"], n_cells=71689, donor_id_available="yes: patient/sample metadata", sex_metadata_available="yes", raw_counts_available="yes: SRA raw reads and processed unnormalized counts", processed_matrix_available="yes: h5ad", cell_type_annotations_available="yes: h5ad", age_available="yes", eligible_for_pseudobulk_sex_DE="yes for control subset only (3M/4F), low-powered", eligible_for_cell_composition_analysis="yes for control subset only (3M/4F), low-powered", analysis_subset="controls: 7 donors (3M/4F); aneurysm: 6 donors (3M/3F)", preferred_minimum_met="yes for control subset", low_power_flag="yes: only 3 male controls", notes="GEO reports 71,689 nuclei across 13 donors. The published control subset contains 39,346 nuclei (secondary report); healthy/control inference must use the 7 non-aneurysm controls only; do not mix aneurysm with controls for healthy sex DE.", source_url=SOURCE_URLS["GSE207784"]),
        row(dataset_id="GSE189795", dataset_alias="normal and dissected ascending aorta scRNA", tissue="ascending aorta", healthy_or_control="mixed: 4 normal controls and 5 acute dissection samples", technology="scRNA-seq", n_donors_total=c189["n"], n_male=c189["male"], n_female=c189["female"], n_cells="not reported in GEO series metadata", donor_id_available="yes: sample-level metadata", sex_metadata_available="yes", raw_counts_available="SRA raw reads; normalized matrix supplementary", processed_matrix_available="yes: normalized counts", cell_type_annotations_available="yes", age_available="yes", eligible_for_pseudobulk_sex_DE="no for control subset: 0 female controls", eligible_for_cell_composition_analysis="no for control subset: 0 female controls", analysis_subset="controls: 4 donors (4M/0F); dissection: 5 donors (5M/0F)", preferred_minimum_met="no", low_power_flag="yes", notes="All 9 GEO samples are male; the 4 healthy controls are all male. This is not eligible for donor-sex comparison.", source_url=SOURCE_URLS["GSE189795"]),
        row(dataset_id="GSE165824", dataset_alias="Deep learning thoracic aorta / SCP1265", tissue="ascending and descending aorta", healthy_or_control="normal rapid-autopsy human aorta", technology="10x snRNA-seq 3' v3", n_donors_total=h165["donors"], n_male=h165["male"], n_female=h165["female"], n_cells=h165["n_cells"], donor_id_available="yes: donor_id in processed h5ad", sex_metadata_available="yes: donor-level sex in processed h5ad obs", raw_counts_available="yes: SRA raw reads", processed_matrix_available="yes: GEO h5ad and SCP1265", cell_type_annotations_available="yes: processed h5ad", age_available="no in GEO/h5ad obs", eligible_for_pseudobulk_sex_DE="no: only 1 male / 2 female donors", eligible_for_cell_composition_analysis="no: below preferred donor minimum", analysis_subset="3 paired normal donors (Ao4 female, Ao8 female, Ao12 male), 54,092 nuclei", preferred_minimum_met="no", low_power_flag="yes", notes="SCP1265 and GSE165824 are the same study, not independent datasets. GEO sample metadata lacks sex, but the deposited h5ad obs contains donor-level sex; this audit used obs metadata only and did not read expression values.", source_url=SOURCE_URLS["GSE165824"]),
    ]

    fields = list(rows[0].keys())
    with (OUT / "01_single_cell_rescue_dataset_audit.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    priority = [
        row(priority_rank=1, dataset_id="GSE202109", recommended_use="primary kidney sex-DE and composition", usable_subset="all 19 healthy living donors", male_donors=9, female_donors=10, cells=27677, preferred_criterion="PASS", priority_reason="Only purpose-built healthy sex-balanced kidney cohort in the audit; donor sex and age available."),
        row(priority_rank=2, dataset_id="GSE183276", recommended_use="secondary kidney scRNA validation", usable_subset="Normal Reference subset", male_donors=7, female_donors=13, cells=21650, preferred_criterion="PASS", priority_reason="Healthy subset has 20 donors and both sexes; retain disease atlas samples separately."),
        row(priority_rank=3, dataset_id="GSE183277", recommended_use="secondary kidney snRNA validation", usable_subset="Normal Reference subset", male_donors=8, female_donors=10, cells=88460, preferred_criterion="PASS", priority_reason="Large healthy reference subset with balanced donor sex; retain disease atlas samples separately."),
        row(priority_rank=4, dataset_id="GSE207784", recommended_use="exploratory vascular backup", usable_subset="7 non-aneurysm controls", male_donors=3, female_donors=4, cells=39346, preferred_criterion="PASS but low power", priority_reason="Control subset meets >=3 per sex but has only 3 male donors; do not mix aneurysm samples."),
        row(priority_rank=5, dataset_id="GSE182416", recommended_use="not eligible for sex-DE; descriptive thyroid cell-type reference only", usable_subset="all 7 samples", male_donors=0, female_donors=7, cells=54726, preferred_criterion="FAIL", priority_reason="All samples female and pathology is PTC/follicular adenoma despite normal-adjacent thyroid label."),
        row(priority_rank=6, dataset_id="GSE165824", recommended_use="descriptive vascular cell-type reference; not powered for sex-DE", usable_subset="3 paired normal donors", male_donors=1, female_donors=2, cells=54092, preferred_criterion="FAIL", priority_reason="Processed h5ad contains donor sex, but only 1 male and 2 female donors; below the preferred minimum."),
        row(priority_rank=7, dataset_id="GSE189795", recommended_use="not eligible for sex-DE", usable_subset="4 normal controls", male_donors=4, female_donors=0, cells="not reported", preferred_criterion="FAIL", priority_reason="All healthy controls are male."),
    ]
    with (OUT / "02_recommended_dataset_priority.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(priority[0].keys()))
        writer.writeheader()
        writer.writerows(priority)

    report = f"""# M5b — URXP02 single-cell rescue audit

## Scope

This audit only verifies dataset availability, donor sex metadata, healthy/control composition, processed/raw data access, annotations, and donor-level eligibility. No candidate-gene expression, differential expression, module analysis, clustering, reannotation, figures, or new epidemiologic analysis was run.

## Decision summary

- **Primary recommendation: GSE202109 kidney.** It contains 19 healthy living kidney donors (9 male, 10 female) and 27,677 cells, satisfying the prespecified minimum of at least three donors per sex.
- **Secondary kidney validation: GSE183276 and GSE183277.** These are mixed kidney atlases, but their deposited `Normal Reference` subsets are sex-balanced (scCv3: 20 donors, 7M/13F, 21,650 cells; snCv3: 18 donors, 8M/10F, 88,460 cells). Disease samples must remain separate.
- **Vascular backup: GSE207784 controls.** The series has 13 donors and 71,689 nuclei overall; the seven non-aneurysm controls are 3M/4F. This meets the minimum numerically but is low-powered.
- **Thyroid rescue failed for sex comparison in GSE182416.** All seven samples are female, and the sample characteristics contain PTC or follicular adenoma pathology. It can provide a cell-type reference, not donor-sex DE.
- **GSE189795 is not eligible for sex comparison.** All four normal controls are male (and all nine series samples are male).
- **SCP1265/GSE165824 is a separate thoracic-aorta study, not a duplicate of GSE207784.** It has three paired normal donors and 54,092 nuclei. The deposited h5ad `obs` metadata identifies two female donors (Ao4, Ao8) and one male donor (Ao12), but this is below the preferred donor threshold and is descriptive only.

## Eligibility rule

The preferred criterion is at least 3 male and 3 female donors in the eligible healthy/control subset. A dataset can be retained below that threshold for descriptive or exploratory purposes, but it is flagged low-powered and is not treated as a confirmatory donor-sex resource.

## Audit provenance

The audit used official GEO family SOFT records for GSE182416, GSE202109, GSE207784, GSE189795, and GSE165824, deposited kidney atlas metadata for GSE183276/GSE183277, and metadata-only inspection of the deposited GSE165824 h5ad. The kidney metadata parser explicitly accounts for the deposited extra cell-barcode field before the published header. The incomplete local GSE207784 h5ad download was not used to infer counts; the 71,689-nucleus total is taken from the official GEO series record. For GSE165824, only h5ad dimensions and donor/cell-type/disease metadata were read; expression values were not accessed.

Generated: {datetime.now(timezone.utc).isoformat()}
"""
    (OUT / "M5B_SINGLE_CELL_RESCUE_AUDIT_REPORT.md").write_text(report, encoding="utf-8")

    input_files = [
        CACHE / "GSE182416_family.soft.gz", CACHE / "GSE202109_family.soft.gz",
        CACHE / "GSE207784_family.soft.gz", CACHE / "GSE189795_family.soft.gz",
        CACHE / "GSE165824_family.soft.gz", CACHE / "GSE165824_ascending_descending_human_aorta_v1.h5ad", CACHE / "GSE183276_metadata.txt.gz",
        CACHE / "GSE183277_metadata.txt.gz",
    ]
    manifest = {
        "analysis": "M5b_single_cell_rescue_audit",
        "scope": "dataset audit only; no expression or candidate-gene analysis",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "preferred_minimum": {"male_donors": 3, "female_donors": 3},
        "source_urls": SOURCE_URLS,
        "input_files_sha256": {str(p.relative_to(REPO)): sha256(p) for p in input_files if p.exists()},
        "output_files": [
            "01_single_cell_rescue_dataset_audit.csv",
            "02_recommended_dataset_priority.csv",
            "M5B_SINGLE_CELL_RESCUE_AUDIT_REPORT.md",
            "manifest.json",
        ],
        "notes": [
            "GSE165824 is the GEO accession for SCP1265 Deep learning enables genetic analysis of the human thoracic aorta.",
            "GSE207784 and GSE165824 are distinct studies.",
            "No donor sex was inferred when the source metadata did not provide it.",
        ],
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
