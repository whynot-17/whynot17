from __future__ import annotations

import re
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SOURCE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(r"C:\Users\21634\Desktop\whynot.docx")
DEST_DIR = ROOT / "outputs" / "manuscript" / "final_submission"
DEST_DIR.mkdir(parents=True, exist_ok=True)
DEST = DEST_DIR / "MCOP_CRC_manuscript_final.docx"

if not SOURCE.exists():
    raise FileNotFoundError(SOURCE)

replacements = [
    (
        "Full gate definitions and candidate-level audit fields are provided in the Supplementary Methods and Supplementary Table X.",
        "Full gate definitions and candidate-level audit fields are provided in the Supplementary Methods and Supplementary Table S8.",
    ),
    (
        "full stage-specific counts are shown in Fig. 2 and Supplementary Table X.",
        "full stage-specific counts are shown in Fig. 2 and Supplementary Tables S1 and S8.",
    ),
    (
        "The complete rubric and signal-level audit are provided in Supplementary Table X.",
        "The complete rubric and signal-level audit are provided in Supplementary Table S4.",
    ),
    (
        "Outcome-blinded actionability: 267 candidates \u2192 87 human-testable chemical-biomarker mappings \u2192 15 unique NHANES biomarker tests.",
        "Outcome-blinded actionability defines the human screening universe. (A) Prespecified gates reduced 267 core environmental chemicals to 87 human-testable chemical–biomarker mappings. CRC outcome statistics remained behind the outcome firewall until the biomarker-test universe was frozen; 27 mappings met the strict D2/C2/T2 rule. (B) The 87 eligible mappings represented 15 unique NHANES biomarker tests, which formed the denominator for BH-FDR correction. Point area reflects the number of eligible chemical mappings and color denotes biological matrix. (C) MiNP, parent DINP and urinary MCOP retained distinct roles. MiNP was a molecular nominee but failed the direct-detectability gate, parent DINP was not a significant Phase 1 hit, and MCOP entered the human screen as a measurable biomarker for a DINP-related exposure axis. Biomarker translation does not imply chemical equivalence or a direct MCOP molecular hit.",
    ),
]

with zipfile.ZipFile(SOURCE, "r") as zin:
    members = {name: zin.read(name) for name in zin.namelist()}

xml_name = "word/document.xml"
xml = members[xml_name].decode("utf-8")

# Word may split visible text across identical runs. The source was generated
# without tracked changes, so merge adjacent text boundaries for the four
# exact, auditable replacements while preserving paragraph/run properties.
xml_visible = xml
for old, new in replacements:
    if old not in xml_visible:
        # Permit XML tags between words by locating the paragraph that contains
        # all distinctive tokens, then collapse its text only through a targeted regex.
        old_tokens = [re.escape(x) for x in old.split()]
        pattern = r"(?:<[^>]+>)*\s*".join(old_tokens)
        m = re.search(pattern, xml_visible)
        if m:
            # Replacement inside arbitrary XML is unsafe; require pre-merged source.
            raise RuntimeError(f"Phrase is fragmented across XML runs and must be merged first: {old[:60]}")
        raise RuntimeError(f"Expected phrase not found: {old}")
    xml_visible = xml_visible.replace(old, new, 1)

members[xml_name] = xml_visible.encode("utf-8")

tmp = DEST.with_suffix(".tmp.docx")
with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_DEFLATED) as zout:
    for name, data in members.items():
        zout.writestr(name, data)
tmp.replace(DEST)

# Content-level submission checks.
with zipfile.ZipFile(DEST, "r") as z:
    final_xml = z.read(xml_name).decode("utf-8")
assert "Supplementary Table X" not in final_xml
assert "Supplementary Table S8" in final_xml
assert "Supplementary Tables S1 and S8" in final_xml
assert "Supplementary Table S4" in final_xml
assert "Outcome-blinded actionability defines the human screening universe" in final_xml
assert "direct MCOP molecular hit" in final_xml

print(DEST)
