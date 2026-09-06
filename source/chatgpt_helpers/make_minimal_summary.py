from pathlib import Path
import json
out=Path(r"E:\chatgpt\qc_dinp_minimal_0_50ns")
a=json.loads((out/'status_minimal_0_50ns.json').read_text(encoding='utf-8'))
b=json.loads((out/'core_internal_and_pocket_residue_summary.json').read_text(encoding='utf-8'))
summary={
 'trajectory':'E:/chatgpt/qc_50ns/production.dcd (snapshot through 50.23 ns; analysis window 0.01-50.00 ns)',
 'definitions':a['definitions'],
 'reference_pocket_residues':b['reference_pocket_residues'],
 'reference_pocket_atoms':b['reference_pocket_atoms'],
 'reference_core_pocket_pairs':b['reference_core_pocket_pairs'],
 'overall_0_50ns':{
   'core_internal_rmsd_A':b['overall_0_50ns']['core_internal_rmsd_mean_A'],
   'core_internal_rmsd_sd_A':b['overall_0_50ns']['core_internal_rmsd_sd_A'],
   'core_internal_rmsd_max_A':b['overall_0_50ns']['core_internal_rmsd_max_A'],
   'core_in_pocket_rmsd_A':a['overall_0_50ns']['core_rmsd_mean_A'],
   'core_in_pocket_rmsd_sd_A':a['overall_0_50ns']['core_rmsd_sd_A'],
   'core_in_pocket_rmsd_max_A':a['overall_0_50ns']['core_rmsd_max_A'],
   'branch1_rmsd_A':a['overall_0_50ns']['branch1_rmsd_mean_A'],
   'branch2_rmsd_A':a['overall_0_50ns']['branch2_rmsd_mean_A'],
   'com_displacement_A':a['overall_0_50ns']['com_displacement_mean_A'],
   'com_displacement_max_A':a['overall_0_50ns']['com_displacement_max_A'],
   'pair_retention':b['overall_0_50ns']['pair_retention_mean'],
   'residue_contact_fraction':b['overall_0_50ns']['residue_contact_fraction_mean'],
   'nearest_core_pocket_distance_A':b['overall_0_50ns']['nearest_distance_mean_A'],
   'nearest_distance_le_4p5_fraction':b['overall_0_50ns']['nearest_le_4p5_fraction'],
 },
 'segment_40_50ns':{
   'core_internal_rmsd_A':b['segment_40_50ns']['core_internal_rmsd_mean_A'],
   'core_in_pocket_rmsd_A':a['overall_0_50ns']['core_rmsd_mean_A'],
   'branch1_rmsd_A':a['overall_0_50ns']['branch1_rmsd_mean_A'],
   'branch2_rmsd_A':a['overall_0_50ns']['branch2_rmsd_mean_A'],
   'com_displacement_A':a['overall_0_50ns']['com_displacement_mean_A'],
   'pair_retention':b['segment_40_50ns']['pair_retention_mean'],
   'residue_contact_fraction':b['segment_40_50ns']['residue_contact_fraction_mean'],
   'nearest_core_pocket_distance_A':b['segment_40_50ns']['nearest_distance_mean_A'],
   'nearest_distance_le_4p5_fraction':b['segment_40_50ns']['nearest_le_4p5_fraction'],
 },
 'exact_50ns':{'core_in_pocket_rmsd_A':a['nearest_50ns']['core_rmsd_A'],'core_internal_rmsd_A':b['nearest_50ns']['core_internal_rmsd_A'],'branch1_rmsd_A':a['nearest_50ns']['branch1_rmsd_A'],'branch2_rmsd_A':a['nearest_50ns']['branch2_rmsd_A'],'com_displacement_A':a['nearest_50ns']['com_displacement_A'],'pair_retention':b['nearest_50ns']['pair_retention'],'residue_contact_fraction':b['nearest_50ns']['residue_contact_fraction'],'nearest_core_pocket_distance_A':b['nearest_50ns']['nearest_distance_A']},
 'interpretation':{
   'ring_integrity':'Aromatic ring is rigid (internal RMSD mean ~0.046 A, max ~0.100 A).',
   'branch_behavior':'Both ester/alkyl arms show large core-aligned RMSD, consistent with flexible torsional rearrangement.',
   'pocket_behavior':'The core moved away from the minimized-reference pose (protein-aligned core RMSD ~5.2 A after 10 ns), while the nearest distance to the reference pocket stayed ~3.5 A and <=4.5 A in ~99.8% of frames. Reference atom-pair retention is low (~11% overall, ~8% at 40-50 ns), so the original contact pattern is not retained even though the ligand remains near the pocket region.',
   'decision':'No evidence for sustained complete dissociation in 0-50 ns, but this is not a clean low-core-RMSD/high-original-contact-retention case. Treat it as an early pose relocation to a neighboring pocket microstate plus flexible arm motion. Inspect the relocated core pose/contact map before committing to the remaining 100 ns production.'
 },
 'files':{'summary_json':'E:/chatgpt/qc_dinp_minimal_0_50ns/minimal_test_summary.json','timeseries_csv':'E:/chatgpt/qc_dinp_minimal_0_50ns/dinp_minimal_timeseries.csv','core_contacts_csv':'E:/chatgpt/qc_dinp_minimal_0_50ns/core_internal_and_pocket_residue_timeseries.csv'}
}
(out/'minimal_test_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8'); print(json.dumps(summary,ensure_ascii=False))
