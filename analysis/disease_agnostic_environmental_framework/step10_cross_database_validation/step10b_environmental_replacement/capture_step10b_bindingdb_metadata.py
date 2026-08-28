#!/usr/bin/env python3
"""Capture BindingDB's official download-page metadata without querying data."""
from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import requests

OUT = Path(__file__).resolve().parent
URL = "https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp"


def main() -> None:
    retrieved = datetime.now(timezone.utc).isoformat()
    meta: dict[str, object] = {"url": URL, "retrieved_utc": retrieved, "status": "unavailable"}
    try:
        response = requests.get(URL, headers={"User-Agent": "whynot17-step10b-environmental/1.0"}, timeout=(10, 45))
        content = response.content
        meta.update({"status_code": response.status_code, "content_length": len(content), "response_sha256": hashlib.sha256(content).hexdigest(), "content_type": response.headers.get("content-type")})
        if response.status_code == 200:
            html_path = OUT / "bindingdb_download_page_snapshot.html"
            html_path.write_bytes(content)
            text = response.text
            meta["release_date_candidates"] = sorted(set(re.findall(r"20\\d{2}[-/]\\d{2}[-/]\\d{2}|20\\d{2}[-/]\\d{2}|updated[^<]{0,80}", text, flags=re.I)))[:30]
            meta["status"] = "captured"
    except Exception as exc:
        meta["error"] = f"{type(exc).__name__}: {exc}"
    (OUT / "STEP10B_E_BINDINGDB_METADATA.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    snapshot_path = OUT / "STEP10B_E_SOURCE_SNAPSHOT.json"
    if snapshot_path.exists():
        snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))
        snapshot.setdefault("sources", {}).setdefault("E2", {})["download_page_snapshot"] = meta
        snapshot["retrieval_finished_utc"] = retrieved
        snapshot_path.write_text(json.dumps(snapshot, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(meta, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
