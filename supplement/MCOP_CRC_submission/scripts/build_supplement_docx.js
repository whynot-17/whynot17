const fs = require('fs');
const path = require('path');
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  ImageRun, PageBreak, Table, TableRow, TableCell, WidthType, BorderStyle,
  ShadingType, Header, Footer, PageNumber, TableOfContents
} = require('docx');

const root = path.resolve(__dirname, '..', '..', '..');
const pkg = path.join(root, 'supplement', 'MCOP_CRC_submission');
const figDir = path.join(pkg, 'figures');
const outDocx = path.join(pkg, 'MCOP_CRC_Supplementary_Information.docx');

const BLUE = '2166AC';
const TEAL = '287D8E';
const RED = 'B2182B';
const DARK = '222222';
const GREY = '5F6B72';
const LIGHT = 'EEF3F5';

function run(text, opts={}) {
  return new TextRun({ text, font: 'Arial', size: opts.size || 20, bold: !!opts.bold,
    italics: !!opts.italics, color: opts.color || DARK, break: opts.break || 0,
    superScript: !!opts.superScript });
}

function para(text, opts={}) {
  const children = Array.isArray(text) ? text : [run(text, opts)];
  return new Paragraph({ children, alignment: opts.align || AlignmentType.JUSTIFIED,
    spacing: { after: opts.after === undefined ? 120 : opts.after, line: opts.line || 276 },
    keepNext: !!opts.keepNext, pageBreakBefore: !!opts.pageBreakBefore });
}

function heading(text, level=HeadingLevel.HEADING_1) {
  return new Paragraph({ text, heading: level, spacing: { before: level === HeadingLevel.HEADING_1 ? 280 : 180, after: 100 },
    keepNext: true });
}

function bullet(text, level=0) {
  return new Paragraph({ children: [run(text)], bullet: { level },
    spacing: { after: 70, line: 250 } });
}

function fig(stem, width=620, height=630) {
  return new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 100, after: 120 },
    children: [new ImageRun({ data: fs.readFileSync(path.join(figDir, stem + '.png')), type: 'png', transformation: { width, height } })] });
}

function cell(text, width, header=false) {
  return new TableCell({ width: { size: width, type: WidthType.DXA },
    shading: header ? { fill: 'DCE8ED', type: ShadingType.CLEAR } : undefined,
    margins: { top: 70, bottom: 70, left: 85, right: 85 },
    children: [new Paragraph({ children: [run(text, { bold: header, size: 18 })], spacing: { after: 0, line: 230 } })] });
}

function simpleTable(headers, rows, widths) {
  const borders = { top:{style:BorderStyle.SINGLE,size:2,color:'AAB6BC'}, bottom:{style:BorderStyle.SINGLE,size:2,color:'AAB6BC'},
    left:{style:BorderStyle.SINGLE,size:2,color:'AAB6BC'}, right:{style:BorderStyle.SINGLE,size:2,color:'AAB6BC'},
    insideHorizontal:{style:BorderStyle.SINGLE,size:1,color:'D6DDE1'}, insideVertical:{style:BorderStyle.SINGLE,size:1,color:'D6DDE1'} };
  return new Table({ width: { size: widths.reduce((a,b)=>a+b,0), type: WidthType.DXA }, borders,
    rows: [new TableRow({ tableHeader:true, children:headers.map((x,i)=>cell(x,widths[i],true)) }),
      ...rows.map(r=>new TableRow({ children:r.map((x,i)=>cell(String(x),widths[i],false)) }))] });
}

const docChildren = [];

docChildren.push(
  new Paragraph({ alignment: AlignmentType.CENTER, spacing:{before:360,after:220}, children:[run('SUPPLEMENTARY INFORMATION',{size:30,bold:true,color:BLUE})] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing:{after:160}, children:[run('A data-first environmental screening framework identifies urinary MCOP as a robust biomarker signal associated with prevalent colorectal cancer',{size:28,bold:true})] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing:{after:280}, children:[run('Submission-ready methods, figures, tables and reproducibility boundary',{size:20,italics:true,color:GREY})] }),
  para([run('Package status: ',{bold:true}),run('frozen to the latest manuscript and validated analysis outputs. ') ,run('Interpretive lock: ',{bold:true,color:RED}),run('cross-sectional association; no causal or mediation claim.',{color:RED})],{align:AlignmentType.CENTER,after:220}),
  simpleTable(['Component','Contents'],[
    ['Supplementary Methods','Outcome firewall, actionability gates, NHANES survey design, robustness audit and transcriptomic analysis'],
    ['Supplementary Figures','Figures S1–S4, each with vector PDF/SVG and 300-dpi PNG'],
    ['Supplementary Tables','Tables S1–S9 in MCOP_CRC_Supplement_Tables.xlsx'],
    ['Source Data','Panel-level data for Figures S1–S4 in MCOP_CRC_Source_Data.xlsx'],
    ['Reproducibility','Independent R/Python survey agreement and a 10 PASS / 1 INFO / 0 FAIL revalidation audit']
  ],[2300,6500]),
  new Paragraph({children:[new PageBreak()]}),
  heading('Contents',HeadingLevel.HEADING_1),
  new TableOfContents('Supplementary Information', { hyperlink: true, headingStyleRange: '1-3' }),
  new Paragraph({children:[new PageBreak()]})
);

docChildren.push(heading('Terminology and interpretation lock'));
docChildren.push(para('The following terminology is used consistently across the main manuscript, Supplementary Information, figures, tables and source data. The lock prevents chemical identity, statistical unit and evidence-level drift during submission formatting.'));
docChildren.push(simpleTable(['Canonical term','Locked meaning / prohibited interpretation'],[
  ['MCOP','Mono(carboxy-isooctyl) phthalate; urinary NHANES biomarker URXCOP used to index a DINP-related exposure axis. It is not represented as a significant direct CTD molecular hit.'],
  ['MiNP','Monoisononyl phthalate; molecular nominee from Phase 1. It remains chemically distinct from MCOP and failed the frozen direct-detectability gate for primary human screening.'],
  ['DINP-related exposure axis','Exposure-axis language linking parent DINP-related chemistry to urinary MCOP biomonitoring. It does not establish parent-compound-specific causality.'],
  ['Unique biomarker test','The statistical multiplicity unit in the unified NHANES screen. Eighty-seven eligible chemical–biomarker mappings represented 15 unique tests.'],
  ['Prevalent CRC','Current or previous colorectal cancer status in cross-sectional NHANES; not incident risk.'],
  ['PPAR/NR remodeling','Paired CRC epithelial disease-state difference. It is not evidence that DINP or MCOP caused the state.'],
  ['LOCO','Leave-one-cycle-out pooled re-estimation. The seven overlapping estimates are not independent replications.']
],[2450,6350]));

docChildren.push(heading('Supplementary Methods'));
docChildren.push(heading('Study architecture and outcome firewall',HeadingLevel.HEADING_2));
docChildren.push(para('The study was organized as four analytically separated layers: hypothesis-agnostic molecular nomination, outcome-blinded human actionability filtering, systematic NHANES biomarker screening, and independent CRC disease-state transcriptomics. CRC odds ratios, confidence intervals and P values were inaccessible to the actionability rules used to define the human-testable universe. Outcome statistics were introduced only after chemical–biomarker mappings had been resolved and collapsed to 15 unique biomarker tests. This outcome firewall was retained in the 267-row matrix through the field PRIORITIZATION_OUTCOME_BLINDED.'));

docChildren.push(heading('Molecular nomination',HeadingLevel.HEADING_2));
docChildren.push(para('Human direct chemical–gene interaction records from the Comparative Toxicogenomics Database were deduplicated at the chemical–gene level and intersected with colorectal-cancer genes from a GeneCards Disorders-scoped search. The principal background universe comprised genes associated with the final core environmental-chemical set; an all-CTD background was used as a sensitivity universe. For each chemical, enrichment was evaluated with the one-sided hypergeometric/Fisher framework and Benjamini–Hochberg false-discovery-rate control. GeneCards relevance information was retained as a weighted-overlap measure. Degree-matched permutation compared each chemical with chemicals of similar CTD gene-set size to reduce research-degree bias. Chemical-interacting gene set is used instead of target set because CTD interactions include expression, activity and other regulatory relationships and do not necessarily indicate direct molecular binding.'));
docChildren.push(para('MiNP was nominated at this molecular stage (rank 24; BH-FDR 0.00346; degree-matched empirical FDR 0.0356). MCOP had 19 human CTD-interacting genes and two CRC overlaps but did not meet the molecular significance/stability threshold (BH-FDR approximately 0.28). Parent DINP likewise was not a significant Phase 1 hit. Accordingly, the molecular layer nominates a DINP/MiNP-related axis; it does not constitute a direct MCOP molecular discovery.'));

docChildren.push(heading('Outcome-blinded actionability gates',HeadingLevel.HEADING_2));
docChildren.push(para('All 267 core environmental chemicals entered a prespecified actionability audit. Sequential gates were applied to each row without CRC outcome information. E required a resolved chemical entity; X required an interpretable human exposure construct; B required a specific or documented human biomarker mapping; D required adequate detectability; C required sufficient multi-cycle coverage; and T required an analyzable NHANES exposure–outcome–covariate frame with survey design variables. The permissive and moderate rules required E=X=B=1, D at least 1, C at least 1 or 2, respectively, and T at least 1. The strict rule required E=X=B=1 and D=C=T=2.'));
docChildren.push(para('The sequential counts were 267 starting chemicals, 259 E-valid, 135 E+X interpretable, 134 with a biomarker, 127 detectable, 124 with adequate coverage, and 87 human-testable mappings; 27 were strict-eligible. Shared biomarkers were collapsed only for statistical testing, not interpreted as chemical equivalence. The 87 mappings corresponded to 15 unique NHANES biomarker tests.'));

docChildren.push(heading('DINP-axis biomarker translation',HeadingLevel.HEADING_2));
docChildren.push(para('MiNP, DINP parent compounds and MCOP were handled as distinct chemical entities. Direct urinary MiNP detectability was 40.7% in the updated actionability audit and failed the primary detectability gate. MCOP was detectable in 98.8% of measured records, covered seven cycles and entered the systematic screen as the urinary biomarker for a DINP-related exposure axis. This translation reflects biomonitoring utility and does not assert molecular identity, metabolite equivalence or parent-specific causal attribution.'));

docChildren.push(heading('NHANES population, outcome and exposure',HeadingLevel.HEADING_2));
docChildren.push(para('The MCOP analysis pooled seven 2-year NHANES cycles from 2005–2006 through 2017–2018. The harmonized frame contained 17,382 records; 12,127 had MCOP measurements. Before complete-case restriction, the frozen CRC-versus-cancer-free frame contained 11,086 participants, including 76 CRC cases and 11,010 controls. The primary complete-case sample contained 9,936 participants, 70 CRC cases and 9,866 controls. CRC was defined from the frozen current/previous cancer questionnaire algorithm. Controls were cancer-free under the primary definition. Urinary MCOP was log2 transformed so the coefficient represents an exposure doubling. Values below the analytical detection limit were handled using the frozen NHANES laboratory-value construction and were not redefined during the Supplement audit.'));

docChildren.push(heading('Complex-survey design and primary model',HeadingLevel.HEADING_2));
docChildren.push(para('Cycle-specific phthalate subsample weights were selected from each NHANES laboratory component and divided by seven. Strata and primary sampling unit identifiers were prefixed by cycle to prevent accidental cross-cycle merging of reused numeric codes. The complete-case design contained 214 PSUs in 105 strata and no singleton strata; the design degrees of freedom were 109. The primary model was CRC status regressed on log2(MCOP), age, sex, race/ethnicity, body-mass index, smoking, poverty-to-income ratio and log2 urinary creatinine. The frozen primary inference used R survey::svyglm with design-based standard errors and design-degrees-of-freedom P values.'));
docChildren.push(para('A separate Python Newton–IRLS estimator with a Taylor PSU sandwich reproduced the coefficient and standard error. Weight division by 10 in a historical construction and the corrected division by seven differed only by a common scaling factor; because the Python fitter normalized weights by their mean, coefficient and variance estimates were numerically unchanged. The corrected seven-cycle weight is used in all submission outputs.'));

docChildren.push(heading('Systematic 15-test screen and multiplicity',HeadingLevel.HEADING_2));
docChildren.push(para('Each of the 15 frozen biomarker tests was analyzed with the same covariate structure, using the analyte-specific survey weight and the cycles in which the biomarker was measured. Urinary biomarkers additionally included log2 urinary creatinine. The primary family comprised exactly 15 tests and was adjusted with the Benjamini–Hochberg procedure. No subset-based re-adjustment was used to promote selected results. MCOP and PFHS were the two FDR-supported tests and both proceeded to the same robustness scorecard.'));

docChildren.push(heading('Uniform robustness scorecard',HeadingLevel.HEADING_2));
docChildren.push(para('Eight prespecified tags summarized multiplicity and stability. F encoded primary support (2 for BH-FDR<0.05, 1 for nominal P<0.05, 0 otherwise). L encoded LOCO stability (2 for same direction with all LOCO confidence intervals excluding 1, 1 for same direction with some intervals crossing 1, 0 for directional instability). C encoded cycle-direction consistency (2 for at least 80%, 1 for 60–79%, 0 for less than 60%). H encoded exposure-by-cycle heterogeneity (2 for interaction P at least 0.10, 1 for 0.05–<0.10, 0 for <0.05). D and T encoded preservation of direction under diagnosis-timing and upper-tail exclusions, with score 2 requiring maximum absolute log-OR change no greater than 0.25. A encoded algorithmic behavior (2 no warning, 1 localized warning with estimable fits, 0 persistent warning or failure). E encoded event information (2 for at least 60 CRC cases, 1 for 30–59 and 0 for fewer than 30). Robust Tier A required F2, L at least 1, C at least 1, D at least 1, T at least 1 and A at least 1; H was retained as a penalty/evidence tag rather than a deletion gate.'));

docChildren.push(heading('MCOP sensitivity and heterogeneity analyses',HeadingLevel.HEADING_2));
docChildren.push(para('Prespecified analyses included age at least 40 years; sex-specific effects and a formal MCOP-by-sex interaction; exclusion of CRC diagnoses less than 1, 2 or 5 years before examination; exclusion of the top 1% and 2.5% of MCOP; creatinine-normalized MCOP; pairwise adjustment for MEHHP, MEOHP, MECPP and MBzP; adjustment for a non-MCOP phthalate burden; seven LOCO pooled models; seven cycle-specific models; and a global MCOP-by-cycle interaction. LOCO models assessed dependence on any single cycle, whereas the interaction assessed equality of effects across cycles.'));
docChildren.push(para('Quartiles were calculated with both unweighted and survey-weighted cut points. Restricted cubic splines used four knots at the 5th, 35th, 65th and 95th percentiles, with survey-weighted knot construction as a sensitivity analysis. The continuous log2 model remained prespecified as primary; categorical and spline analyses were secondary exposure-shape characterization.'));

docChildren.push(heading('Assay and calendar-cycle audit',HeadingLevel.HEADING_2));
docChildren.push(para('For every cycle, the audit recorded the NHANES codebook, analytical platform, cycle-specific lower limit of detection, proportion above the detection limit, weighted MCOP distribution, creatinine distribution, CRC event count, case and control MCOP medians, weighted CRC prevalence and demographic composition. MCOP was measured by HPLC-ESI-MS/MS methods across all seven cycles. The 2011–2012 cycle had 100% detectability and an LLOD of 0.2 ng/mL; therefore, its discordant effect estimate could not be attributed to gross non-detection. In that cycle the raw complete-case case median was lower than the control median.'));

docChildren.push(heading('CRC transcriptomic disease-state analyses',HeadingLevel.HEADING_2));
docChildren.push(para('Transcriptomic analyses were independent of the NHANES association and did not contain measured DINP or MCOP perturbation. The primary single-cell analysis used the frozen CELLxGENE Census release 2025-11-08 and the verified official dataset-level H5AD source for the eligible paired CRC dataset. All Census queries were restricted to primary data. Individual cells were aggregated to donor-level pseudobulk; cells were not treated as independent statistical replicates. Tumor-derived epithelial was used unless malignant status was directly supported by source annotation.'));
docChildren.push(para('The prespecified PPAR/nuclear-receptor core comprised PPARA, PPARD, PPARG, NR1I2, NR1I3, NR1H2 and NR1H3. RELA and STAT3 formed a separate inflammatory module; the nine-gene composite was not treated as a primary mechanism score because the two modules changed in opposite directions. Paired donor differences were tested with two-sided Wilcoxon signed-rank tests and adjusted across score definitions with BH-FDR. Independent definitions included receptor-only and nuclear-receptor partner modules, KEGG PPAR signaling, Reactome nuclear-receptor and lipid-metabolism pathways, Hallmark metabolic programs, an enterocyte differentiation program, and DoRothEA PPARA/PPARG regulon activities.'));
docChildren.push(para('Within-state analyses required at least 20 cells in both tumor and normal conditions for a donor. Enterocyte-like, secretory-like and other epithelial annotations were analyzed separately. Parallel compartment analyses evaluated epithelial, endothelial, fibroblast and myeloid donor-level contrasts. The evidence-tier lock classified direct paired human disease-state observations as E3, donor-level state associations as E2, external toxicology-supported candidate links as E1, and the untested MCOP/DINP-to-CRC epithelial-state bridge as E0. No formal mediation analysis was claimed, and the causal exposure-to-state arrow was prohibited.'));

docChildren.push(heading('Reproducibility and software',HeadingLevel.HEADING_2));
docChildren.push(para('The frozen analysis was revalidated in an isolated environment. The final audit comprised 11 checks: 10 PASS, 1 INFO and 0 FAIL. The INFO item records that a live staged Census cache completed 35 of 36 donors; the final paired epithelial analysis instead used the separately validated official H5AD and was unaffected. The independent R survey implementation used survey version 4.5. The local CELLxGENE environment used cellxgene-census 1.17.0 and tiledbsoma 1.17.1. Raw NHANES XPT files, large H5AD files and raw-expression caches are not duplicated in the Git repository; their provenance and hashes are retained in the reproducibility inventory.'));

docChildren.push(heading('Supplementary Results notes'));
docChildren.push(heading('The actionability framework fixed the multiplicity denominator before outcome analysis',HeadingLevel.HEADING_2));
docChildren.push(para('The largest loss occurred at exposure interpretation: 124 candidates first failed X, compared with 8 at E, 1 at B, 7 at D, 3 at C and 37 at T. Eighty-seven mappings remained eligible and represented 15 unique biomarker tests. The DINP-axis transition was explicit: MiNP retained the strongest relevant molecular nomination but failed direct detectability, whereas MCOP entered as the measurable urinary biomarker.'));
docChildren.push(heading('MCOP was the sole Robust Tier A signal',HeadingLevel.HEADING_2));
docChildren.push(para('Two tests passed the 15-test BH-FDR threshold: PFHS (OR 0.624, 95% CI 0.471–0.828; BH-FDR 0.0219) and MCOP (OR 1.246, 95% CI 1.077–1.440; BH-FDR 0.0248). MCOP had fingerprint F2|L2|C2|H0|D2|T2|A1|E2 and was the only Robust Tier A signal. PFHS was Tier B because the analysis contained 32 CRC cases and retained a persistent algorithmic warning. The competing PFHS result is reported rather than suppressed.'));
docChildren.push(heading('The pooled MCOP estimate was reproducible but temporally heterogeneous',HeadingLevel.HEADING_2));
docChildren.push(para('R survey::svyglm produced OR 1.245507 (95% CI 1.077309–1.439966; design-df P=0.003311). The independent Python implementation produced the same OR with absolute log-OR difference 1.07×10−13; direction and confidence-interval conclusion agreed. All seven LOCO estimates were positive and excluded the null, but the global MCOP-by-cycle interaction was significant (F-test P=0.00598). Six of seven cycle-specific point estimates exceeded 1; 2011–2012 was the discordant cycle.'));
docChildren.push(heading('CRC transcriptomics supported state-specific remodeling, not an exposure mechanism',HeadingLevel.HEADING_2));
docChildren.push(para('The frozen seven-gene PPAR/NR score was lower in paired tumor-derived epithelium (36 donors; median Δ −0.419; BH-FDR 9.30×10−7). Twelve of 13 estimable alternative definitions decreased; PPARD regulon activity was not estimable. The pattern was compartment-specific: epithelial PPAR/NR decreased, endothelial and fibroblast contrasts were near null, and myeloid PPAR/NR increased. Within epithelial annotations, enterocyte-like cells decreased (24 paired donors; median Δ −0.191; BH-FDR 5.26×10−4), whereas secretory-like cells increased modestly (27 donors; median Δ +0.068; BH-FDR 0.0181). These results support epithelial-state-specific remodeling and leave the environmental-to-state link untested.'));

docChildren.push(new Paragraph({children:[new PageBreak()]}));
docChildren.push(heading('Supplementary Figures'));

docChildren.push(heading('Supplementary Figure S1. Auditable actionability filtering and multiplicity denominator.',HeadingLevel.HEADING_2));
docChildren.push(para('A, Sequential outcome-blinded filtering of 267 core environmental chemicals to 87 human-testable chemical–biomarker mappings. The blue firewall denotes that CRC outcome statistics were not used before the 15-test universe was frozen. The strict rule retained 27 mappings. B, Number of candidates first failing each gate; a candidate is counted at its earliest failure. C, The 87 eligible mappings represented 15 unique NHANES biomarker tests; point area reflects the number of eligible chemical mappings and color denotes biological matrix. D, Locked chemical-identity interpretation for MiNP, parent DINP and MCOP. Biomarker translation does not imply chemical equivalence.'));
docChildren.push(fig('Figure_S1_actionability_audit',620,703));

docChildren.push(heading('Supplementary Figure S2. Full human screen and uniform robustness audit.',HeadingLevel.HEADING_2));
docChildren.push(para('A, Effect landscape for all 15 frozen biomarker tests. The x axis is log2 OR per exposure doubling, the y axis is −log10 BH-FDR, point area scales with CRC events and pale horizontal intervals show 95% confidence intervals. MCOP and PFHS passed the 15-test BH-FDR threshold; URXMOH was nominally positive but did not pass FDR. B, Prespecified robustness fingerprints. F, multiplicity support; L, LOCO stability; C, cycle-direction consistency; H, exposure-by-cycle homogeneity; D, diagnosis-timing stability; T, upper-tail stability; A, algorithmic behavior; E, event information. C, Independent R and Python implementations of the primary MCOP model.'));
docChildren.push(fig('Figure_S2_human_screen_robustness',620,636));

docChildren.push(heading('Supplementary Figure S3. Calendar-cycle exposure and assay audit.',HeadingLevel.HEADING_2));
docChildren.push(para('A, Survey-weighted MCOP distribution by NHANES cycle. The shaded region denotes Q1–Q3; solid and dotted series denote the median and P95. B, Percentage above the analytical detection limit and codebook LLOD for each cycle. C, Cycle-specific ORs per MCOP doubling with 95% confidence intervals, CRC event counts and raw complete-case case/control MCOP median ratios. The pooled OR and confidence band are shown for context. The 2011–2012 estimate was discordant despite 100% detectability; the global MCOP-by-cycle interaction P value was 0.006. Cycle-specific analyses are underpowered and are not interpreted as seven independent replications.'));
docChildren.push(fig('Figure_S3_cycle_exposure_audit',620,624));

docChildren.push(heading('Supplementary Figure S4. Definition robustness, state specificity and causal boundary.',HeadingLevel.HEADING_2));
docChildren.push(para('A, Paired tumor-minus-normal median differences across independent PPAR/nuclear-receptor, metabolic and regulon score definitions; displayed values are BH-FDR. Twelve of 13 estimable definitions decreased. B, Frozen PPAR/NR contrasts within epithelial annotations and across major cellular compartments. Point area scales with paired donor count. C, Evidence-tier lock. E3 denotes direct paired human disease-state observations, E2 donor-level state associations, E1 an external-toxicology-supported candidate link and E0 the untested MCOP/DINP-to-CRC epithelial-state bridge. The E0 causal arrow is prohibited.'));
docChildren.push(fig('Figure_S4_ppar_state_evidence',620,578));

docChildren.push(new Paragraph({children:[new PageBreak()]}));
docChildren.push(heading('Supplementary Tables'));
docChildren.push(para('All tables are supplied in MCOP_CRC_Supplement_Tables.xlsx. Each worksheet is also available as a CSV in the reproducibility package. The workbook is intended to be the authoritative tabular supplement; panel-level figure source data are supplied separately in MCOP_CRC_Source_Data.xlsx.'));
docChildren.push(simpleTable(['Table','Title and content'],[
  ['Table S1','Sample selection and actionability attrition: frozen NHANES counts and 267→87→15 filtering summary.'],
  ['Table S2','Primary MCOP complex-survey model, weight-scaling audit, age-restricted model and independent R/Python implementation comparison.'],
  ['Table S3','Complete 15-axis human biomarker screen with analytic N, CRC events, OR, 95% CI, P value, BH-FDR and screen rank.'],
  ['Table S4','Uniform eight-domain robustness scorecard for all 15 biomarker tests, including warnings and final robustness tier.'],
  ['Table S5','MCOP sensitivity analyses: LOCO, age, sex, diagnosis timing, upper-tail exclusion, creatinine normalization and pairwise co-exposure models.'],
  ['Table S6','Cycle-specific estimates, global exposure-by-cycle interaction and restricted cubic spline tests.'],
  ['Table S7','Cycle-level MCOP distribution, detectability, LLOD, assay platform and case/control exposure summary.'],
  ['Table S8','Complete 267-row actionability matrix with entity, mapping, detectability, coverage, testability, eligibility and outcome-firewall fields.'],
  ['Table S9','Transcriptomic evidence lock, DINP-axis candidate bridge audit, PPAR/NR definitions, within-state contrasts and causal-status fields.']
],[1650,7150]));

docChildren.push(heading('Data and code availability'));
docChildren.push(para('The repository contains analysis scripts, frozen small outputs, figure code, supplementary tables, source-data workbooks and audit records. Large raw NHANES downloads, official H5AD files and Census/raw-expression caches are intentionally excluded from version control. The submission package records their provenance and local hash inventory. The repository should be archived to a persistent release before acceptance, with the release identifier inserted into the final Data Availability statement.'));

docChildren.push(heading('Submission consistency checklist'));
[
  'Main and Supplement use 267 starting chemicals, 87 eligible mappings and 15 unique biomarker tests.',
  'MCOP primary complete-case N=9,936 with 70 CRC cases.',
  'Primary R survey OR=1.245507; 95% CI 1.077309–1.439966; design-df P=0.003311.',
  'The unified screen reports both FDR-supported tests: MCOP and PFHS.',
  'MCOP is the sole Robust Tier A biomarker; significant cycle heterogeneity remains visible.',
  'MiNP, DINP and MCOP are not treated as chemically equivalent.',
  'PPAR/NR results use donor-level inference; cells are not independent replicates.',
  'No direct exposure-to-CRC causal or epithelial-state mediation claim is made.',
  'The live 35/36-donor staged Census cache is an INFO item only; the final analysis uses the validated official H5AD.',
  'Supplementary figures have PDF, SVG and 300-dpi PNG exports and matching source-data sheets.'
].forEach(x=>docChildren.push(bullet(x)));

const doc = new Document({
  creator: 'MCOP–CRC analysis team', title: 'MCOP–CRC Supplementary Information',
  description: 'Submission-ready supplementary methods, figures, tables and reproducibility lock',
  styles: {
    default: { document: { run: { font: 'Arial', size: 20, color: DARK }, paragraph: { spacing: { line: 276 } } } },
    paragraphStyles: [
      { id:'Title', name:'Title', basedOn:'Normal', next:'Normal', quickFormat:true, run:{font:'Arial',size:34,bold:true,color:BLUE}, paragraph:{alignment:AlignmentType.CENTER,spacing:{after:220}} },
      { id:'Heading1', name:'Heading 1', basedOn:'Normal', next:'Normal', quickFormat:true, run:{font:'Arial',size:28,bold:true,color:BLUE}, paragraph:{spacing:{before:280,after:120},keepNext:true,outlineLevel:0} },
      { id:'Heading2', name:'Heading 2', basedOn:'Normal', next:'Normal', quickFormat:true, run:{font:'Arial',size:23,bold:true,color:TEAL}, paragraph:{spacing:{before:200,after:90},keepNext:true,outlineLevel:1} },
      { id:'Heading3', name:'Heading 3', basedOn:'Normal', next:'Normal', quickFormat:true, run:{font:'Arial',size:20,bold:true,color:DARK}, paragraph:{spacing:{before:150,after:70},keepNext:true,outlineLevel:2} }
    ]
  },
  sections: [{
    properties: { page: { size: { width: 11906, height: 16838 }, margin: { top: 1080, right: 960, bottom: 1080, left: 960 } } },
    headers: { default: new Header({ children:[new Paragraph({ alignment:AlignmentType.RIGHT, children:[run('MCOP–CRC | Supplementary Information',{size:16,color:GREY})] })] }) },
    footers: { default: new Footer({ children:[new Paragraph({ alignment:AlignmentType.CENTER, children:[run('Page ',{size:16,color:GREY}), new TextRun({children:[PageNumber.CURRENT],font:'Arial',size:16,color:GREY})] })] }) },
    children: docChildren
  }]
});

Packer.toBuffer(doc).then(buffer => {
  fs.mkdirSync(pkg, { recursive: true });
  fs.writeFileSync(outDocx, buffer);
  console.log(outDocx);
});
