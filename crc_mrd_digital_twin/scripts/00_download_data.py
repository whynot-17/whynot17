"""Download the complete public supplementary package for Chen et al. (2021).

Source: Springer Nature Figshare, article DOI 10.1186/s13045-021-01089-z.
The script downloads all associated Additional files and verifies each Figshare MD5.
"""
from __future__ import annotations
import csv, hashlib, json, urllib.request
from pathlib import Path

ARTICLE_DOI = "10.1186/s13045-021-01089-z"
SEARCH_URL = "https://api.figshare.com/v2/articles/search"
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data_raw"


def get_json(url: str):
    req = urllib.request.Request(url, headers={"User-Agent": "crc-mrd-audit/1.0"})
    with urllib.request.urlopen(req, timeout=60) as response:
        return json.load(response)


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    payload = json.dumps({"resource_doi": ARTICLE_DOI, "limit": 100}).encode()
    req = urllib.request.Request(
        SEARCH_URL, data=payload,
        headers={"Content-Type": "application/json", "User-Agent": "crc-mrd-audit/1.0"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=60) as response:
        items = json.load(response)
    items = [x for x in items if x.get("resource_doi") == ARTICLE_DOI]
    expected = {f"Additional file {i} " for i in range(1, 19)}
    found = {x["title"][:len(prefix)] for x in items for prefix in expected if x["title"].startswith(prefix)}
    if len(items) != 18 or len(found) != 18:
        raise RuntimeError(f"Expected Additional files 1-18; found {len(items)} items")

    manifest = []
    for item in sorted(items, key=lambda x: int(x["title"].split("Additional file ")[1].split()[0])):
        meta = get_json(item["url"])
        for file in meta.get("files", []):
            target = RAW_DIR / file["name"]
            expected_md5 = file.get("supplied_md5", "").lower()
            digest = hashlib.md5()
            if target.exists():
                with target.open("rb") as existing:
                    while chunk := existing.read(1024 * 1024):
                        digest.update(chunk)
            actual = digest.hexdigest()
            if actual != expected_md5:
                request = urllib.request.Request(file["download_url"], headers={"User-Agent": "crc-mrd-audit/1.0"})
                digest = hashlib.md5()
                with urllib.request.urlopen(request, timeout=120) as response, target.open("wb") as out:
                    while chunk := response.read(1024 * 1024):
                        out.write(chunk)
                        digest.update(chunk)
                actual = digest.hexdigest()
            if actual != expected_md5:
                raise RuntimeError(f"MD5 mismatch for {file['name']}: {actual} != {expected_md5}")
            manifest.append({
                "additional_file": int(meta["title"].split("Additional file ")[1].split()[0]),
                "title": meta["title"], "figshare_id": meta["id"], "doi": meta["doi"],
                "source_doi": meta["resource_doi"], "filename": file["name"],
                "bytes": file["size"], "md5": actual, "download_url": file["download_url"],
                "license": meta.get("license", {}).get("name"),
            })
    with (RAW_DIR / "figshare_manifest.json").open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
    with (RAW_DIR / "figshare_manifest.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=manifest[0].keys())
        writer.writeheader()
        writer.writerows(manifest)
    print(f"Downloaded and MD5-verified {len(manifest)} files from {len(items)} Figshare records.")


if __name__ == "__main__":
    main()
