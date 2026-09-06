import sqlite3
from pathlib import Path

p = Path("work/geo/annotation/hugene10/hugene10sttranscriptcluster.db/inst/extdata/hugene10sttranscriptcluster.sqlite")
ids = {
    "BBOX1": "8424", "SLC22A5": "6584", "CPT1A": "1374", "CPT1B": "1375", "CPT2": "1376",
    "SLC25A20": "788", "ACADM": "34", "ACADVL": "37", "HADHA": "3030", "HADHB": "3032",
    "ECHS1": "1892", "ETFA": "2108", "ETFB": "2109", "ETFDH": "2110", "PPARGC1A": "10891",
    "KLF5": "688", "FABP6": "2172",
}
con = sqlite3.connect(p)
for symbol, gene_id in ids.items():
    probes = [row[0] for row in con.execute("select probe_id from probes where gene_id=?", (gene_id,))]
    print(symbol, gene_id, len(probes), probes[:8])
con.close()
