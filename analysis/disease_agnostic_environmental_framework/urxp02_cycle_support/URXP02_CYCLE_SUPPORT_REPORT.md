# URXP02 cycle support audit

This is a read-only alignment check using the frozen URXP02 primary/robustness cycle metadata and the frozen outcome inventory. No association model was refit, no robustness estimate was recomputed, and no external source was accessed.

## Exposure cycle support

URXP02 (`2-hydroxynaphthalene` in the frozen exposure registry) is listed in exactly 8 exposure cycles: 2001-2002, 2003-2004, 2005-2006, 2007-2008, 2009-2010, 2011-2012, 2013-2014, 2015-2016. However, the outcome inventory marks thyroid disease as unavailable in 2001-2002. Hypertension is available there. The thyroid `LOCO:2001-2002` refit is a no-op (its N and interaction estimate equal the full robustness fit), confirming that no thyroid outcome records from that cycle entered the fitted sample. The true common outcome-definable support is therefore 7 cycles: 2003-2004, 2005-2006, 2007-2008, 2009-2010, 2011-2012, 2013-2014, 2015-2016.

## Outcome definitions

- Thyroid disease: Age >=20. Case: MCQ160M=Yes. Control: MCQ160M=No.
- Hypertension: Age >=20. Case: BPQ020=Yes. Control: BPQ020=No.

## Sex-specific cycle counts

The inventory records below are descriptive cycle-level counts. The original frozen availability rule was based on pooled sex-specific case counts and at least four definable cycles; no new per-cycle case threshold was invented here.

- 2001-2002: thyroid disease male=0, female=0; hypertension male=727, female=911; both outcome records available=False; both-sex cases present=False; both endpoints contribute observations=False.
- 2003-2004: thyroid disease male=95, female=403; hypertension male=828, female=921; both outcome records available=True; both-sex cases present=True; both endpoints contribute observations=True.
- 2005-2006: thyroid disease male=86, female=373; hypertension male=760, female=814; both outcome records available=True; both-sex cases present=True; both endpoints contribute observations=True.
- 2007-2008: thyroid disease male=102, female=452; hypertension male=996, female=1128; both outcome records available=True; both-sex cases present=True; both endpoints contribute observations=True.
- 2009-2010: thyroid disease male=120, female=488; hypertension male=1048, female=1144; both outcome records available=True; both-sex cases present=True; both endpoints contribute observations=True.
- 2011-2012: thyroid disease male=118, female=393; hypertension male=973, female=1030; both outcome records available=True; both-sex cases present=True; both endpoints contribute observations=True.
- 2013-2014: thyroid disease male=97, female=504; hypertension male=988, female=1157; both outcome records available=True; both-sex cases present=True; both endpoints contribute observations=True.
- 2015-2016: thyroid disease male=151, female=471; hypertension male=1001, female=1069; both outcome records available=True; both-sex cases present=True; both endpoints contribute observations=True.

The 8-cycle simultaneous-definition check therefore fails at 2001–2002 because thyroid disease is unavailable. The corrected 7-cycle common-support check passes: both outcomes are available, both sexes have cases, and both outcomes contribute observations in each of the 7 common cycles. This does not by itself establish causal validity or fixed-FDR confirmation; it sets the correct cycle scope for any follow-up claim.
