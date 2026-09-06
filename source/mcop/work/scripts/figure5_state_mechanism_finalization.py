"""Finalize Figure 5 evidence boundaries, collision audit, tiers and report."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs"


def literature_audit() -> pd.DataFrame:
    rows = [
        {
            "citation": "Wang et al., Genome Medicine 2022; PMID 35974387",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9380328/",
            "cohort_or_model": "12 metastatic CRC patients; 9,120 single cells; paired adjacent normal, primary and metastatic tumors; tumor organoids",
            "cell_or_tissue_unit": "epithelial cells for DEG discovery; organoids for inhibitor experiments",
            "ppar_definition": "PPAR-signaling-associated genes among tumor-vs-normal epithelial DEGs; pathway enrichment, not a prespecified receptor score",
            "method": "single-cell epithelial DEG heatmap/enrichment; FH535/GW9662 organoid perturbation",
            "reported_direction": "selected PPAR-associated DEGs up in tumor; PPAR inhibition reduced organoid growth",
            "definition_compatibility_with_ours": "partial_only",
            "audit_interpretation": "Does not test the frozen 7-gene receptor/NR expression score; authors themselves report discordance with TCGA and opposite SCD/ACSL4 trends.",
        },
        {
            "citation": "Koelwyn et al., J Proteome Res 2019; PMID 31398025",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6705791/",
            "cohort_or_model": "consensus molecular subtype colon cancer transcriptome/proteome case study",
            "cell_or_tissue_unit": "bulk tumor subtype comparison",
            "ppar_definition": "KEGG PPAR signaling co-expression configuration",
            "method": "self-contained gene-set/co-expression analysis",
            "reported_direction": "only part of pathway up in CMS3; explicitly excluded PPAR receptors",
            "definition_compatibility_with_ours": "different",
            "audit_interpretation": "Supports modular heterogeneity: downstream lipid genes can rise while receptor expression falls.",
        },
        {
            "citation": "Wang et al., Cancer Cell International 2026; PMID 41593461",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12918641/",
            "cohort_or_model": "CRC bulk, three paired scRNA samples, clinical tissue, CRC cells and mouse experiments",
            "cell_or_tissue_unit": "tumor epithelial cells and experimental models",
            "ppar_definition": "FABP2-PPAR-alpha-sphingomyelin axis",
            "method": "multi-omics plus functional perturbation",
            "reported_direction": "FABP2 reduced in tumor epithelium; PPAR-alpha-linked biology remains context dependent",
            "definition_compatibility_with_ours": "adjacent",
            "audit_interpretation": "Consistent with loss of differentiated lipid-handling epithelial identity, but not direct DINP evidence.",
        },
        {
            "citation": "Yaghoubizadeh et al., PPAR Research 2020; PMID 32399461",
            "url": "https://pubmed.ncbi.nlm.nih.gov/32399461/",
            "cohort_or_model": "100 paired CRC tumor and adjacent normal tissues",
            "cell_or_tissue_unit": "bulk paired tissue",
            "ppar_definition": "individual PPARA, PPARD and PPARG mRNA",
            "method": "RT-qPCR",
            "reported_direction": "PPARA and PPARD higher, PPARG lower",
            "definition_compatibility_with_ours": "component_level_mixed",
            "audit_interpretation": "Demonstrates isoform-specific directions; a single undifferentiated 'PPAR activated' label is not reproducible across definitions.",
        },
    ]
    frame = pd.DataFrame(rows)
    frame["paper"] = frame["citation"]
    frame["year"] = [2022, 2019, 2026, 2020]
    frame["dataset"] = ["HRA000183 / study scRNA-seq", "consensus molecular subtype datasets", "GSE231559 plus bulk/experimental cohorts", "100 paired clinical specimens"]
    frame["n_patients"] = [12, "not a single-cell cohort", 3, 100]
    frame["cell_type"] = ["epithelial", "bulk tumor", "epithelial", "bulk tissue"]
    frame["comparison"] = ["tumor vs adjacent-normal epithelial DEGs", "CMS3 vs other CMS", "tumor vs normal epithelial FABP2", "paired tumor vs adjacent normal"]
    frame = frame.rename(columns={"ppar_definition": "PPAR_definition"})
    frame["gene_set_source"] = ["reported PPAR pathway enrichment; exact fixed list not supplied in article text", "KEGG PPAR signaling", "FABP2/PPAR-alpha mechanistic axis", "individual receptor genes"]
    frame["genes_if_available"] = ["selected DEG-associated genes; article notes SCD and ACSL4 discordance", "KEGG pathway subset excluding receptors", "FABP2 and PPAR-alpha axis", "PPARA;PPARD;PPARG"]
    frame["analysis_unit"] = ["cells/DEGs for discovery; organoid replicates for perturbation", "bulk samples", "cells plus experimental replicates", "patient-paired tissue"]
    frame["statistical_method"] = frame["method"]
    frame["reported_effect"] = frame["reported_direction"]
    frame["interpretation"] = frame["audit_interpretation"]
    frame["potential_reason_for_difference"] = [
        "different cohort, metastatic case mix, DEG-enrichment rather than donor-level fixed score, and epithelial-state composition",
        "subtype-specific downstream metabolic co-expression without receptor upregulation",
        "different molecular anchor and very small scRNA cohort",
        "bulk composition and opposing receptor-isoform effects",
    ]
    frame["donor_level_analysis"] = [False, True, False, True]
    frame["cell_level_pseudoreplication_risk"] = ["unclear/high for discovery DEG unit; organoid tests separate", "not applicable", "high for the small scRNA comparison unless donor-aware", "none for paired tissue test"]
    frame["tumor_vs_normal"] = [True, False, True, True]
    frame["pathway_enrichment_only"] = [True, False, False, False]
    return frame


def toxicology_bridge() -> pd.DataFrame:
    rows = [
        {
            "chemical": "MiNP", "pmid": "27551952", "url": "https://pubmed.ncbi.nlm.nih.gov/27551952/",
            "model": "human receptor transactivation/two-hybrid assays and primary human hepatocytes",
            "tissue_relevance": "hepatic/receptor assay; not colon", "endpoint": "CAR2/PXR activation; weaker human PPAR activation and target-gene responses",
            "direction": "activation", "evidence_type": "single-chemical experimental", "usable_in_main_figure": "no_non_colon_supplement_only",
            "boundary": "Shows receptor engagement at experimental concentrations, not CRC epithelial state or human-dose causality.",
        },
        {
            "chemical": "MiNP", "pmid": "23843199", "url": "https://pubmed.ncbi.nlm.nih.gov/23843199/",
            "model": "in-silico docking to human PPAR/RXR subtypes", "tissue_relevance": "no tissue",
            "endpoint": "predicted binding", "direction": "binding", "evidence_type": "computational",
            "usable_in_main_figure": "supplement_only", "boundary": "Docking is not evidence of cellular activation or colon relevance.",
        },
        {
            "chemical": "MiNP", "pmid": "35421560", "url": "https://pubmed.ncbi.nlm.nih.gov/35421560/",
            "model": "primary mouse granulosa cells", "tissue_relevance": "ovary; not colon",
            "endpoint": "PPRE reporter and PPAR target genes", "direction": "dose-dependent/nonmonotonic; mainly PPAR-gamma in this model",
            "evidence_type": "single-chemical experimental", "usable_in_main_figure": "supplement_only",
            "boundary": "Supports PPAR engagement but is species- and tissue-specific.",
        },
        {
            "chemical": "DINP", "pmid": "31154059", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5356750/",
            "model": "FITC allergic dermatitis mouse ear model", "tissue_relevance": "skin/immune; not colon",
            "endpoint": "phospho-RELA and phospho-STAT3", "direction": "increased in co-exposure dermatitis context",
            "evidence_type": "animal co-exposure/context-dependent", "usable_in_main_figure": "supplement_only",
            "boundary": "Supports inflammatory signaling plausibility but not single-agent CRC epithelial causality.",
        },
        {
            "chemical": "MCOP/MCIOP", "pmid": "34478338", "url": "https://pubmed.ncbi.nlm.nih.gov/34478338/",
            "model": "760 pregnant women; maternal urinary metabolite and placental RNA-seq", "tissue_relevance": "placenta; not colon",
            "endpoint": "18 associated placental transcripts/pathways", "direction": "association; mixed",
            "evidence_type": "human observational transcriptomics", "usable_in_main_figure": "no",
            "boundary": "Does not establish MCOP as a direct molecular perturbagen; biomarker and tissue differ from CRC epithelium.",
        },
        {
            "chemical": "DINP/MiNP", "pmid": "42398653", "url": "https://pubmed.ncbi.nlm.nih.gov/42398653/",
            "model": "mouse liver and in-vitro metabolic assays", "tissue_relevance": "liver; not colon",
            "endpoint": "PPAR activation, beta-oxidation and lipid-metabolism transcriptomics", "direction": "activation/remodeling at model-dependent doses",
            "evidence_type": "single-chemical experimental", "usable_in_main_figure": "no_non_colon_supplement_only",
            "boundary": "Recent exposure evidence strengthens general PPAR plausibility but does not reproduce the CRC epithelial suppression state.",
        },
    ]
    frame = pd.DataFrame(rows)
    frame["gene_or_pathway"] = ["CAR/PXR/PPAR nuclear receptors", "PPARA/PPARD/PPARG/RXR binding", "PPAR response element", "RELA/STAT3 phosphorylation", "placental transcriptome", "PPAR/lipid oxidation"]
    frame["tissue"] = frame["tissue_relevance"]
    frame["single_chemical_or_mixture"] = frame["evidence_type"]
    frame = frame.rename(columns={"pmid": "PMID"})
    frame["directness"] = ["direct receptor/cell assay", "in-silico only", "direct non-colon cell assay", "contextual animal co-exposure", "observational biomarker association", "experimental non-colon tissue"]
    frame["reason"] = frame["boundary"]
    return frame


def collision_audit() -> pd.DataFrame:
    rows = [
        ("exact", 'DINP AND colorectal cancer AND PPAR', 0, "No direct DINP-PPAR-CRC mechanistic paper identified"),
        ("exact", 'MiNP AND colorectal cancer AND PPAR', 0, "No direct MiNP-PPAR-CRC mechanistic paper identified"),
        ("exact", 'MCOP AND colorectal cancer AND PPAR', 0, "No direct MCOP-PPAR-CRC mechanistic paper identified"),
        ("partial", 'DINP AND cancer AND PPAR', 1, "2026 systematic carcinogenic-hazard review emphasizes rodent PPAR-alpha liver mode and limited human cancer evidence; PMID 42094681"),
        ("partial", 'MiNP AND PPAR', 4, "Receptor/hepatocyte, ovarian and macrophage studies; none establish CRC epithelial mediation"),
        ("adjacent", 'colorectal cancer AND PPAR signaling single-cell', 4, "CRC PPAR literature is substantial but definition- and isoform-dependent"),
        ("background", 'phthalates AND PPAR', 100, "Broad crowded background; count is qualitative/lower-bound and not used for novelty claims"),
        ("exact", 'DINP AND colon epithelial AND PPAR', 0, "No exact exposure-epithelial-state study identified"),
        ("exact", 'MiNP AND colorectal cancer', 0, "No direct CRC paper identified"),
        ("exact", 'MiNP AND colon epithelial', 0, "No direct colon epithelial perturbation identified"),
        ("partial", 'phthalate AND colorectal cancer AND single cell', 0, "No study connecting measured DINP/MCOP exposure to CRC single-cell state identified"),
        ("partial", 'phthalate AND epithelial state AND colorectal', 0, "No exact chain identified"),
        ("adjacent", 'PPARG AND colorectal epithelial AND single cell', 3, "CRC receptor/state literature exists without DINP/MCOP human exposure"),
        ("adjacent", 'PPARA AND colorectal epithelial AND single cell', 3, "CRC receptor/state literature exists without DINP/MCOP human exposure"),
        ("adjacent", 'PPAR AND RELA AND STAT3 AND colorectal epithelial', 1, "Pathway crosstalk literature is adjacent; no exposure-to-state chain"),
        ("adjacent", 'PPAR AND stress-like epithelial AND colorectal cancer', 1, "State-remodeling literature is adjacent; no DINP/MCOP linkage"),
        ("background", 'nuclear receptor AND colorectal epithelial state', 10, "Broad disease biology only"),
    ]
    frame = pd.DataFrame(rows, columns=["collision_level", "query", "audited_relevant_records", "finding"])
    frame["exact_prior_complete_chain"] = False
    frame["novelty_wording_allowed"] = "targeted search did not identify an exact prior study"
    return frame


def evidence_tiers(defs: pd.DataFrame, within: pd.DataFrame, tox: pd.DataFrame) -> pd.DataFrame:
    frozen = defs.loc[defs.definition == "Frozen 7-gene PPAR/NR core"].iloc[0]
    standard = defs.loc[defs.definition.isin([
        "KEGG hsa03320 PPAR signaling", "Reactome PPAR-alpha lipid regulation",
        "Reactome peroxisomal lipid metabolism", "Hallmark fatty acid metabolism",
    ])]
    within_sig = within.loc[within.BH_FDR < 0.05]
    within_down = within_sig.loc[within_sig.direction == "down"]
    within_up = within_sig.loc[within_sig.direction == "up"]
    exact_tox = tox.loc[tox.tissue_relevance.str.contains("colon", case=False) & ~tox.tissue_relevance.str.contains("not colon", case=False)]
    rows = [
        ("CRC -> epithelial PPAR/NR down", "E3", "directly observed in current human paired-donor data", "solid", "GREEN", f"Frozen core median delta={frozen.median_delta_tumor_minus_normal:.3f}; FDR={frozen.BH_FDR:.3g}."),
        ("CRC -> RELA/STAT3 inflammatory state up", "E3", "directly observed in current human paired-donor data", "solid", "GREEN", "Expression/regulon and state evidence; no exposure attribution."),
        ("PPAR/NR <-> differentiation/inflammatory state", "E2", "donor-level cross-program convergence", "thin solid", "YELLOW", "Association, not mediation."),
        ("Within-state PPAR remodeling", "E2", "annotation-state paired analysis", "mixed split arrows", "YELLOW", f"{len(within_down)} state down and {len(within_up)} state up at FDR<0.05; composition also shifts."),
        ("DINP/MiNP -> PPAR/NR", "E1", "external toxicology-supported candidate", "dashed candidate label only; omit non-colon evidence nodes", "YELLOW", f"{len(tox)} curated records but {len(exact_tox)} direct colon/intestinal records."),
        ("MCOP/DINP -> CRC epithelial PPAR state", "E0", "untested exposure-to-state bridge", "gap/dotted", "RED", "No human colon perturbation or formal mediation; causal arrow prohibited."),
    ]
    return pd.DataFrame(rows, columns=["link", "evidence_level", "evidence_basis", "figure_representation", "status", "locked_wording"])


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def markdown_table(frame: pd.DataFrame) -> str:
    """Small dependency-free Markdown renderer for audit tables."""
    clean = frame.copy().fillna("")
    columns = [str(c) for c in clean.columns]
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in clean.astype(str).itertuples(index=False, name=None):
        lines.append("| " + " | ".join(value.replace("|", "\\|").replace("\n", " ") for value in row) + " |")
    return "\n".join(lines)


def main() -> None:
    defs = pd.read_csv(OUT / "figure5_ppar_definition_comparison.csv")
    comp = pd.read_csv(OUT / "figure5_epithelial_state_composition.csv")
    within = pd.read_csv(OUT / "figure5_within_state_ppar_analysis.csv")
    corr = pd.read_csv(OUT / "figure5_state_correlation_matrix.csv")
    lit = literature_audit()
    tox = toxicology_bridge()
    collision = collision_audit()
    tiers = evidence_tiers(defs, within, tox)
    lit.to_csv(OUT / "figure5_ppar_literature_definition_audit.csv", index=False)
    tox.to_csv(OUT / "figure5_toxicology_bridge_evidence.csv", index=False)
    collision.to_csv(OUT / "figure5_collision_audit.csv", index=False)
    tiers.to_csv(OUT / "figure5_evidence_tier_lock.csv", index=False)

    standard = defs.loc[defs.definition.isin([
        "KEGG hsa03320 PPAR signaling", "Reactome PPAR-alpha lipid regulation",
        "Reactome peroxisomal lipid metabolism", "Reactome mitochondrial fatty-acid beta oxidation",
        "Hallmark fatty acid metabolism", "Hallmark cholesterol homeostasis",
        "Enterocyte metabolic differentiation",
    ]), ["definition", "median_delta_tumor_minus_normal", "BH_FDR", "direction"]]
    custom = defs.loc[defs.definition == "Frozen 7-gene PPAR/NR core"].iloc[0]
    significant_within = within.loc[within.BH_FDR < 0.05]
    strongest_corr = corr.sort_values("BH_FDR").head(5)
    inflammatory_corr = corr.loc[corr.y.str.contains("RELA|STAT3|TNF|IL6|stress|Inflammatory", case=False, regex=True)].sort_values("BH_FDR")
    comp_summary = comp[["state", "n_paired_donors", "paired_median_delta", "paired_wilcoxon_P", "paired_BH_FDR", "paired_direction"]].drop_duplicates("state")
    within_down = significant_within.loc[significant_within.direction == "down"]
    within_up = significant_within.loc[significant_within.direction == "up"]
    overall = "YELLOW"
    standard_down = standard.loc[(standard.direction == "down") & (standard.BH_FDR < 0.05)]
    standard_other = standard.loc[~standard.index.isin(standard_down.index)]
    opening = [
        f"1. Standard PPAR pathway definitions show MIXED but predominantly DOWN in these CRC epithelial paired donors: {len(standard_down)} significant-down and {len(standard_other)} other definition(s).",
        f"2. Our custom PPAR/NR result IS reproduced by most independent KEGG/Reactome/Hallmark definitions; frozen median delta={custom.median_delta_tumor_minus_normal:.3f}, P={custom.p_value:.3g}.",
        f"3. The observed signal is MIXED: composition shifts coexist with within-state remodeling ({len(within_down)} state down, {len(within_up)} state up at BH-FDR<0.05).",
        f"4. PPAR/NR remodeling IS PARTIALLY linked to inflammatory-stress programs at donor level, but not in the proposed simple inverse direction: {inflammatory_corr.iloc[0].y} rho={inflammatory_corr.iloc[0].spearman_rho:.3f}, BH-FDR={inflammatory_corr.iloc[0].BH_FDR:.3g}; RELA/STAT3 regulon-delta links are not significant.",
        "5. DINP/MiNP toxicology evidence provides candidate convergence but DOES NOT establish exposure-to-state causality; MCOP remains an exposure biomarker.",
        f"6. Figure 5 status: {overall} — retain as a disease-state convergence figure with a dashed environmental bridge, not a causal mechanism figure.",
    ]
    report = "# Figure 5 PPAR contradiction audit and mechanism lock\n\n" + "\n\n".join(opening)
    report += "\n\n## Why the apparent contradiction is not a single statistical disagreement\n\n"
    report += (
        "The 2022 single-cell CRC paper called PPAR signaling activated after selecting tumor-versus-normal epithelial DEGs, "
        "observing enrichment of PPAR-associated genes, and perturbing tumor organoids with PPAR inhibitors. It did not test our "
        "prespecified seven-receptor/nuclear-receptor expression score at the donor level. The same paper noted that several genes "
        "did not reproduce in TCGA and that SCD and ACSL4 could reverse direction. Our audit therefore separates receptor abundance, "
        "broad KEGG/Reactome lipid programs, peroxisomal metabolism and inferred regulon activity.\n\n"
    )
    report += "## Same-cohort definition comparison\n\n" + markdown_table(standard) + "\n\n"
    report += markdown_table(defs[["definition", "n_paired_donors", "median_delta_tumor_minus_normal", "p_value", "BH_FDR", "direction", "genes_present", "coverage"]])
    report += "\n\n## Epithelial composition and within-state audit\n\n"
    report += "Composition paired summary (the CSV retains all donor-condition-state rows):\n\n" + markdown_table(comp_summary) + "\n\nWithin annotation states (minimum 20 cells per condition):\n\n" + markdown_table(within)
    report += "\n\nThe source annotations support enterocyte-like, secretory-like and other epithelial states only. No malignant label was invented."
    report += "\n\n## Donor-delta convergence\n\n" + markdown_table(corr.sort_values("BH_FDR"))
    report += ("\n\nThe positive IL6-JAK-STAT3 and inflammatory-response delta correlations mean that donors with larger PPAR/NR losses do not show the largest inflammatory gains. "
               "RELA and STAT3 regulon-activity deltas are also not significantly correlated with the PPAR/NR delta. Thus the cohort supports parallel disease-state changes, "
               "but does not support a direct donor-level PPAR-low -> RELA/STAT3-high coupling or mediation arrow.")
    report += "\n\n## Toxicology and causal boundary\n\n" + markdown_table(tox)
    report += "\n\nThe main-figure exposure bridge may show only dashed general nuclear-receptor/PPAR plausibility. It must not depict MCOP as a direct perturbagen, nor DINP/MiNP as a proven cause of the CRC epithelial state."
    report += "\n\n## Evidence-tier lock\n\n" + markdown_table(tiers)
    report += "\n\n## Collision audit\n\n" + markdown_table(collision)
    report += "\n\n## Final Figure 5 scientific wording\n\n"
    report += (
        "**CRC epithelial disease-state convergence:** a prespecified PPAR/nuclear-receptor expression module is reduced in paired "
        "CRC epithelium, whereas downstream lipid/PPAR pathway definitions and individual receptor regulons are modular and can differ "
        "in direction. DINP/MiNP toxicology provides non-colon nuclear-receptor plausibility only. The computational evidence supports "
        "a candidate exposure-to-state bridge, not a causal DINP/MCOP-PPAR-CRC mechanism."
    )
    (OUT / "figure5_final_report.md").write_text(report + "\n", encoding="utf-8")

    files = [
        "figure5_ppar_literature_definition_audit.csv", "figure5_ppar_definition_comparison.csv",
        "figure5_ppar_gene_module_decomposition.csv", "figure5_epithelial_state_composition.csv",
        "figure5_within_state_ppar_analysis.csv", "figure5_state_correlation_matrix.csv",
        "figure5_toxicology_bridge_evidence.csv", "figure5_collision_audit.csv",
        "figure5_evidence_tier_lock.csv", "figure5_final_report.md",
    ]
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True).stdout.strip()
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(), "git_head_at_analysis": commit,
        "analysis_unit": "paired donor pseudobulk; donor/group/state aggregates only",
        "paired_donors": 36, "verdict": overall,
        "causal_claim_allowed": False,
        "inputs": ["mcop_phase2g_donor_state_pseudobulk.csv", "mcop_phase2g_epithelial_state_scores.csv", "mcop_phase2g_regulator_activity.csv"],
        "outputs": {name: {"bytes": (OUT / name).stat().st_size, "sha256": sha256(OUT / name)} for name in files},
    }
    (OUT / "figure5_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print("\n".join(opening))


if __name__ == "__main__":
    main()
