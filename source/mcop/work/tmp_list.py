from pathlib import Path
for f in ['work/gene_sets/h.all.v2026.1.Hs.symbols.gmt','work/gene_sets/c2.cp.reactome.v2026.1.Hs.symbols.gmt']:
    print('\n',f)
    for line in Path(f).read_text(encoding='utf-8').splitlines():
        name=line.split('\t',1)[0]
        if any(x in name for x in ['HALLMARK_','FERROPTOSIS','NRF2','GLUTATHIONE','AUTOPHAGY','PURINE','PYRIMIDINE','UNFOLDED_PROTEIN','OXIDATIVE_PHOSPHORYLATION','FATTY_ACID','CHOLESTEROL','ABC_FAMILY','XENOBIOTIC','APOPTOSIS','TGF_BETA','TNFA_SIGNALING','IL6_JAK_STAT3','EMT','GLYCOLYSIS','DNA_REPAIR']):
            print(name)
