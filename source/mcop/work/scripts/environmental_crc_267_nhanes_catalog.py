"""Build a local, reproducible catalog of CDC NHANES laboratory files.

The catalog is the source index for the 267-chemical human actionability
audit. It is intentionally broader than the existing phthalate-only local
data and records every environmental laboratory domain exposed by the CDC
continuous NHANES laboratory pages.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
from bs4 import BeautifulSoup


ROOT = Path(__file__).resolve().parents[2]
CATALOG_DIR = ROOT / "work" / "nhanes_phase2a" / "catalogs"
OUTPUT_DIR = ROOT / "outputs"
CATALOG_OUT = OUTPUT_DIR / "environmental_crc_267_nhanes_lab_catalog.csv"
ENV_OUT = OUTPUT_DIR / "environmental_crc_267_nhanes_environmental_lab_catalog.csv"

CYCLE_MAP = {
    1999: "1999-2000", 2001: "2001-2002", 2003: "2003-2004",
    2005: "2005-2006", 2007: "2007-2008", 2009: "2009-2010",
    2011: "2011-2012", 2013: "2013-2014", 2015: "2015-2016",
    2017: "2017-2018",
}

ENVIRONMENTAL_TERMS = [
    "metal", "perfluoro", "polyfluoro", "pesticide", "phthalate", "plasticizer",
    "pah", "polyaromatic", "polycyclic aromatic", "volatile organic", "voc",
    "dioxin", "furan", "pcb", "flame retard", "phenol", "paraben", "cotinine",
    "organophosphate", "carbamate", "organochlorine", "environmental pollut",
    "air pollut", "water pollut", "brominated flame",
]


def absolute_url(href: str) -> str:
    if href.startswith("http"):
        return href
    return "https://wwwn.cdc.gov" + href


def parse_file(path: Path) -> list[dict[str, str]]:
    match = re.search(r"lab_(\d{4})\.html$", path.name)
    if not match:
        return []
    year = int(match.group(1))
    cycle = CYCLE_MAP[year]
    soup = BeautifulSoup(path.read_text(encoding="utf-8", errors="replace"), "html.parser")
    rows = []
    for tr in soup.select("tr"):
        cells = tr.find_all("td")
        if len(cells) < 3:
            continue
        title = cells[0].get_text(" ", strip=True)
        links = {str(a.get_text(" ", strip=True)).lower(): absolute_url(a.get("href")) for a in tr.find_all("a") if a.get("href")}
        doc = next((url for label, url in links.items() if "doc" in label), "")
        data = next((url for label, url in links.items() if "data" in label and url.lower().endswith(".xpt")), "")
        if not doc and not data:
            continue
        rows.append({
            "cycle": cycle,
            "cycle_begin_year": year,
            "laboratory_title": title,
            "doc_url": doc,
            "data_url": data,
            "data_file": Path(data).name if data else "",
            "environmental_domain": any(term in title.lower() for term in ENVIRONMENTAL_TERMS),
            "catalog_source": str(path),
        })
    return rows


def main() -> None:
    files = sorted(CATALOG_DIR.glob("lab_*.html"))
    if len(files) != len(CYCLE_MAP):
        raise FileNotFoundError(f"Expected {len(CYCLE_MAP)} CDC catalog pages, found {len(files)}")
    rows = [row for path in files for row in parse_file(path)]
    catalog = pd.DataFrame(rows).drop_duplicates(["cycle", "data_url"]).sort_values(["cycle_begin_year", "laboratory_title"])
    env = catalog.loc[catalog["environmental_domain"]].copy()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    catalog.to_csv(CATALOG_OUT, index=False)
    env.to_csv(ENV_OUT, index=False)
    print({"catalog_rows": len(catalog), "environmental_rows": len(env), "unique_environmental_data_files": int(env["data_url"].nunique()), "output": str(CATALOG_OUT)})


if __name__ == "__main__":
    main()
