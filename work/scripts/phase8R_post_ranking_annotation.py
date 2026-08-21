#!/usr/bin/env python
"""Phase 8-R post-ranking annotation.

This script is intentionally separate from the phenotype screen.  The four
drugs below were frozen by the local GDSC phenotype-only gate first; this
script then records a dated, auditable regulatory/novelty review.  These
annotations never change the biological ranking.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs"
AUDIT_DATE = "2026-08-21"

# Sources were reviewed after the biological shortlist was frozen.  The
# targeted audit is not a proof of absence; it records what was found and
# prevents a candidate from being described as novel without checking.
ANNOTATIONS = {
    "SN-38": {
        "aliases": "7-ethyl-10-hydroxycamptothecin",
        "primary_targets": "TOP1-DNA complex",
        "regulatory_status": "not a standalone approved medicine; active metabolite of approved irinotecan",
        "oncology_status": "established oncology cytotoxic metabolite",
        "clinical_usability": "indirect via irinotecan; not a repurposable non-oncology drug",
        "nononcology_eligible": False,
        "novelty_class": "C_direct_overlap",
        "novelty_summary": "SN-38 is an established irinotecan metabolite in CRC; CRC acquired SN-38/oxaliplatin resistance models are already described.",
        "post_ranking_decision": "NO-GO as Drug X; retain as an irinotecan/TOP1 comparator",
        "source_urls": "https://dailymed.nlm.nih.gov/dailymed/fda/fdaDrugXsl.cfm?setid=de634e80-328a-4473-b395-0172683a2411&type=display;https://pubmed.ncbi.nlm.nih.gov/25759163/",
    },
    "AZD2014": {
        "aliases": "vistusertib",
        "primary_targets": "mTORC1/mTORC2",
        "regulatory_status": "investigational; FDA not approved for any disease in the reviewed trial record",
        "oncology_status": "oncology investigational agent",
        "clinical_usability": "clinical-trial exposure exists, but not an approved repurposing drug",
        "nononcology_eligible": False,
        "novelty_class": "C_direct_crc_overlap",
        "novelty_summary": "AZD2014 has already been tested in colon cancer with irinotecan/SN-38 and compared with FOLFOX/FOLFIRI in preclinical studies.",
        "post_ranking_decision": "NO-GO as novel Drug X; retain as an mTOR pharmacologic comparator",
        "source_urls": "https://clinicaltrials.gov/study/NCT03071874;https://pmc.ncbi.nlm.nih.gov/articles/PMC6826690/",
    },
    "JQ1": {
        "aliases": "+JQ1 / BET inhibitor",
        "primary_targets": "BET bromodomains, especially BRD4",
        "regulatory_status": "research chemical/prototype inhibitor; not clinically usable as JQ1",
        "oncology_status": "preclinical oncology probe",
        "clinical_usability": "low; metabolic instability and no approved indication",
        "nononcology_eligible": False,
        "novelty_class": "B_crc_overlap_not_direct_oxa",
        "novelty_summary": "JQ1 has established CRC studies and CRC combinations with irinotecan/topoisomerase-I-directed therapy, but the targeted audit did not identify a direct OXA-R CRC validation for this exact state-first claim.",
        "post_ranking_decision": "NO-GO as repurposable Drug X; retain as a BET mechanism comparator",
        "source_urls": "https://pubmed.ncbi.nlm.nih.gov/34998911/;https://pubmed.ncbi.nlm.nih.gov/32981006/;https://pubmed.ncbi.nlm.nih.gov/38084912/",
    },
    "AZD1332": {
        "aliases": "AZD-1332",
        "primary_targets": "NTRK1/TrkA, NTRK2/TrkB, NTRK3/TrkC",
        "regulatory_status": "research/investigational compound; no approved indication located",
        "oncology_status": "toolbox/investigational kinase inhibitor",
        "clinical_usability": "low; no established clinical dosing or approved indication located",
        "nononcology_eligible": False,
        "novelty_class": "A_adjacent_only",
        "novelty_summary": "The targeted audit found the compound identity/Trk pharmacology and adjacent computational CRC drug-response mentions, but no direct CRC OXA-R mechanistic study for this state-first claim.",
        "post_ranking_decision": "CONDITIONAL mechanistic lead only; NO-GO as approved non-oncology Drug X",
        "source_urls": "https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId=7700;https://pmc.ncbi.nlm.nih.gov/articles/PMC12648302/",
    },
}


def main() -> None:
    shortlist = pd.read_csv(OUT / "phase8R_biological_shortlist.csv")
    cov = pd.read_csv(OUT / "phase8R_covariate_adjusted_drug_associations.csv")
    rows = []
    for rec in shortlist.to_dict("records"):
        drug = rec["DRUG_NAME"]
        ann = ANNOTATIONS[drug]
        c = cov[cov["DRUG_NAME"] == drug]
        rows.append({
            "DRUG_NAME": drug,
            "drug_key": rec["drug_key"],
            "biological_rank": rec["biological_rank"],
            "mean_median_background_rho": rec["mean_median_background_rho"],
            "min_median_background_rho": rec["min_median_background_rho"],
            "GDSC1_global_empirical_q_value": rec["GDSC1_global_empirical_q_value"],
            "GDSC2_global_empirical_q_value": rec["GDSC2_global_empirical_q_value"],
            "primary_targets": ann["primary_targets"],
            "aliases": ann["aliases"],
            "regulatory_status": ann["regulatory_status"],
            "oncology_status": ann["oncology_status"],
            "clinical_usability": ann["clinical_usability"],
            "nononcology_eligible": ann["nononcology_eligible"],
            "covariate_rows": len(c),
            "covariate_median_raw_rho": c["raw_rho"].median() if len(c) else float("nan"),
            "covariate_median_adjusted_partial_rho": c["adjusted_partial_rho"].median() if len(c) else float("nan"),
            "novelty_class": ann["novelty_class"],
            "novelty_summary": ann["novelty_summary"],
            "post_ranking_decision": ann["post_ranking_decision"],
            "annotation_date": AUDIT_DATE,
            "annotation_source_urls": ann["source_urls"],
        })
    annotation = pd.DataFrame(rows).sort_values("biological_rank")
    annotation.to_csv(OUT / "phase8R_drug_regulatory_annotation.csv", index=False)

    nononcology = annotation[annotation["nononcology_eligible"]].copy()
    nononcology["nononcology_shortlist_status"] = "eligible_after_biological_freeze"
    nononcology.to_csv(OUT / "phase8R_nononcology_shortlist.csv", index=False)

    final = annotation.copy()
    final["final_DrugX_eligible"] = final["nononcology_eligible"]
    final["final_selection_status"] = final["final_DrugX_eligible"].map({True: "eligible_pending_further_validation", False: "not_eligible_as_approved_nononcology_DrugX"})
    final.to_csv(OUT / "phase8R_final_drugX_candidates.csv", index=False)

    novelty_lines = [
        "# Phase 8-R novelty audit",
        "",
        f"Audit date: {AUDIT_DATE}. Biological ranking was frozen before this audit.",
        "",
        "Search policy: each shortlisted compound was checked with the standardized concepts `drug + colorectal cancer`, `drug + oxaliplatin`, `drug + oxaliplatin resistance`, and `drug + FOLFOX`. A targeted search is not proof of absence; A/B/C labels describe the evidence found in this audit.",
        "",
        "| Rank | Drug | Class | Audit conclusion |",
        "|---:|---|---|---|",
    ]
    for r in annotation.to_dict("records"):
        novelty_lines.append(f"| {int(r['biological_rank'])} | {r['DRUG_NAME']} | {r['novelty_class']} | {r['novelty_summary']} |")
    novelty_lines += [
        "",
        "A = no direct CRC/OXA-R mechanistic hit located in the targeted audit, but this does not establish absolute novelty; B = CRC overlap without a direct exact OXA-R claim; C = direct CRC/chemotherapy or established metabolite overlap that materially weakens a new Drug X claim.",
        "",
        "## Sources",
        "",
        "- AZD1332 identity and TrkA/B/C pharmacology: https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId=7700",
        "- AZD2014 colon-cancer/irinotecan/SN-38 and FOLFOX/FOLFIRI preclinical study: https://pmc.ncbi.nlm.nih.gov/articles/PMC6826690/",
        "- JQ1 CRC studies and topoisomerase-I-directed combination: https://pubmed.ncbi.nlm.nih.gov/34998911/ and https://pubmed.ncbi.nlm.nih.gov/32981006/",
        "- SN-38 as irinotecan active metabolite and CRC chemoresistance models: https://dailymed.nlm.nih.gov/dailymed/fda/fdaDrugXsl.cfm?setid=de634e80-328a-4473-b395-0172683a2411&type=display and https://pubmed.ncbi.nlm.nih.gov/25759163/",
    ]
    (OUT / "phase8R_novelty_audit.md").write_text("\n".join(novelty_lines) + "\n", encoding="utf-8")

    go_lines = [
        "# Phase 8-R Go/No-Go",
        "",
        "## Decision",
        "",
        "GO for the phenotype-first result as a reproducible computational discovery layer; NO-GO for an approved non-oncology Drug X at this stage.",
        "",
        f"The stringent biological screen produced {len(annotation)} cross-database shortlist drugs. The post-ranking filter produced {len(nononcology)} approved/non-oncology candidates.",
        "",
        "## Candidate disposition",
        "",
        "| Rank | Drug | Decision | Role retained |",
        "|---:|---|---|---|",
    ]
    for r in annotation.to_dict("records"):
        role = "mechanistic comparator" if r["DRUG_NAME"] != "AZD1332" else "conditional mechanistic lead"
        go_lines.append(f"| {int(r['biological_rank'])} | {r['DRUG_NAME']} | {r['post_ranking_decision']} | {role} |")
    go_lines += [
        "",
        "## Interpretation",
        "",
        "The data support a reproducible acquired-OXA-R-like state-conditioned collateral-sensitivity signal, not a universal approved repurposing hit. PRISM and CTRPv2 were unavailable locally and remain missing rather than negative. Before wet-lab prioritization, the next computational gate should be independent pharmacogenomic replication or direct testing of the highest-value mechanistic lead in an external dataset; no candidate should be described as clinically usable from this screen alone.",
    ]
    (OUT / "phase8R_go_nogo.md").write_text("\n".join(go_lines) + "\n", encoding="utf-8")

    report = (OUT / "phase8R_report.md").read_text(encoding="utf-8")
    report += "\n\n## Post-ranking regulatory and novelty audit\n\n"
    report += f"The four-drug biological shortlist was audited on {AUDIT_DATE} after ranking freeze. No candidate met the approved/non-oncology Drug X criterion. SN-38 is retained as an irinotecan/TOP1 comparator; AZD2014 and JQ1 as oncology mechanism comparators; AZD1332 as a conditional Trk mechanistic lead, not as an approved repurposing drug. See `phase8R_drug_regulatory_annotation.csv`, `phase8R_novelty_audit.md` and `phase8R_go_nogo.md`.\n"
    (OUT / "phase8R_report.md").write_text(report, encoding="utf-8")

    manifest_path = OUT / "phase8R_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["status"] = "biological_screen_and_post_ranking_audit_completed_no_approved_nononcology_hit"
    manifest["biological_shortlist"] = annotation["DRUG_NAME"].tolist()
    manifest["approved_nononcology_candidates"] = nononcology["DRUG_NAME"].tolist()
    manifest["post_ranking_audit_date"] = AUDIT_DATE
    manifest["post_ranking_steps_pending"] = ["independent pharmacogenomic replication or direct mechanistic validation"]
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
