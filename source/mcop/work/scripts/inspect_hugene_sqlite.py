import sqlite3
from pathlib import Path

p = Path("work/geo/annotation/hugene10/hugene10sttranscriptcluster.db/inst/extdata/hugene10sttranscriptcluster.sqlite")
con = sqlite3.connect(p)
print([row[0] for row in con.execute("select name from sqlite_master where type='table'")])
for name in ["probes", "gene_info", "chromosome", "probe", "map_counts"]:
    try:
        print(name, con.execute(f"pragma table_info({name})").fetchall()[:12])
    except sqlite3.OperationalError:
        pass
for name in ["accessions", "metadata", "map_metadata"]:
    print(name, con.execute(f"pragma table_info({name})").fetchall()[:20])
print("accessions sample", con.execute("select * from accessions limit 3").fetchall())
print("map names", con.execute("select map_name from map_counts").fetchall())
print("metadata", con.execute("select * from metadata limit 10").fetchall())
print(con.execute("select sql from sqlite_master where type='table' and name in ('probes','accessions')").fetchall())
print("probes sample", con.execute("select * from probes limit 5").fetchall())
print("indexes", con.execute("select name, tbl_name, sql from sqlite_master where type='index'").fetchall()[:10])
print("views", con.execute("select name, sql from sqlite_master where type='view'").fetchall())
con.close()
