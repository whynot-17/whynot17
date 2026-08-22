# Phase 8-M：Meldonium named-candidate challenge

## Prespecification

Meldonium was tested as a named candidate after the Phase 8-R2 signature, thresholds and model universe were frozen. It was not reintroduced into the unbiased drug ranking and no prior candidate list was used.

## Alias-level public-data audit

```text
   platform  metadata_or_response_match  n_matching_records  n_crc_response_records matched_aliases matched_names
       GDSC                       False                   0                       0                              
      PRISM                       False                   0                       0                              
     CTRPv2                       False                   0                       0                              
LINCS/L1000                       False                   0                       0                              
```

Challenge status: `no_direct_response_unavailable_not_negative`.

No Meldonium direct response or LINCS perturbation record was found across the audited local GDSC, PRISM, CTRPv2 and LINCS/L1000 inputs. Therefore no rho, IC50/AUC, percentile, RRS or OXA-R-state drug score was fabricated. The absence is an availability result, not a biological negative.

## Mechanism-consistency audit

```text
             evidence_layer                                                 feature                                          metric                                 value                                                                                                                                        interpretation
OXA-R transcriptional state                                         carnitine_entry                       n_down/n_up; median_delta                        5/1; -0.104503               carnitine-entry is directionally down in 5/6 models; FAO is split 3/3. This is transcript-level state evidence, not flux or dependency.
OXA-R transcriptional state                                       FAO_mitochondrial                       n_down/n_up; median_delta                        3/3; -0.067597               carnitine-entry is directionally down in 5/6 models; FAO is split 3/3. This is transcript-level state evidence, not flux or dependency.
   R3 functional dependency                                                 SLC22A5 r3_rank; median_rho; global_q; leave_HCT116_rho   13863; 0.024459; 0.961814; 0.006421                                       No gene passes the frozen final shortlist; the carnitine/FAO axis lacks stable universal functional dependency.
   R3 functional dependency                                                   CPT1A r3_rank; median_rho; global_q; leave_HCT116_rho 11408; -0.040224; 0.922193; -0.060245                                       No gene passes the frozen final shortlist; the carnitine/FAO axis lacks stable universal functional dependency.
   R3 functional dependency                                                    CPT2 r3_rank; median_rho; global_q; leave_HCT116_rho    5989; 0.081818; 0.751836; 0.173232                                       No gene passes the frozen final shortlist; the carnitine/FAO axis lacks stable universal functional dependency.
   R3 functional dependency                                                   BBOX1 r3_rank; median_rho; global_q; leave_HCT116_rho  10538; 0.047150; 0.904552; -0.013709                                       No gene passes the frozen final shortlist; the carnitine/FAO axis lacks stable universal functional dependency.
     Meldonium pharmacology Meldonium→BBOX1/GBBD inhibition→carnitine availability↓                                  directionality    not a reversal of carnitine-entry↓ The broad direction is same-side rather than opposite-side to the dominant carnitine-entry transcript state; this weakens a universal reversal claim.
Direct public drug response                                               Meldonium             GDSC/PRISM/CTRPv2/LINCS alias audit   0 direct response/signature records                                                   No public pharmacogenomic score or percentile can be calculated; this is unavailable, not negative.
```

## Decision

Meldonium is not rescued by the current dry-lab data. Carnitine-entry is down in 5/6 acquired OXA-R contrasts, but FAO is only 3/3 directional and SLC22A5/CPT1A/CPT2/BBOX1 do not form a stable R3 functional dependency. Its established pharmacology is primarily BBOX/γ-butyrobetaine hydroxylase and OCTN2-related lowering of carnitine availability, rather than a simple direct CPT1A inhibitor ([Dambrova et al., PMID 26850121](https://pubmed.ncbi.nlm.nih.gov/26850121/)). Therefore the broad direction is not a clean reversal of the dominant carnitine-entry state.

The only decisive next test is a matched parental/OXA-R mini-screen. Recommended primary comparison: HCT116-P/OXA-R and DLD1-P/OXA-R or HT29-P/OXA-R, with Meldonium dose-response and pre-specified selectivity index `SI = IC50_parental / IC50_OXA-R`. Suggested interpretation: SI > 1.5-2 in at least two independent backgrounds = revive; SI approximately 1 = deprioritize; SI < 1 = No-Go.

## Source note

Alias registry follows PubChem CID 123868 and NCBI MeSH entry terms: https://pubchem.ncbi.nlm.nih.gov/compound/123868 and https://www.ncbi.nlm.nih.gov/mesh/67050147. Raw response files remain local and are not committed.
