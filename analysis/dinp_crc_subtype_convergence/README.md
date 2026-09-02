# DINP–CRC anatomical subtype convergence

This is a prespecified follow-up of the accessible three-source DINP–CRC
intersection. It asks whether the 81 genes in the current CTD ×
(GeneCards/Open Targets) intersection are preferentially represented in
anatomical CRC concepts in Open Targets.

The analysis uses fixed source-native concepts:

- right-sided anchor: `ascending colon cancer` (`MONDO_0002238`);
- strict left-sided anchor: `sigmoid colon cancer` (`MONDO_0001464`);
- expanded left-sided sensitivity: `rectosigmoid carcinoma` (`MONDO_0002424`),
  combined with sigmoid colon cancer.

Open Targets does not provide a single exact `right-sided colorectal cancer`
versus `left-sided colorectal cancer` concept in this query. Therefore these
outputs are an anatomical knowledge-base localization analysis, not a
patient-level right-versus-left tumor expression analysis and not an
epidemiologic effect modification test.

Run from the repository root:

```powershell
C:\Users\21634\anaconda3\python.exe analysis/dinp_crc_subtype_convergence/run_dinp_crc_subtype_convergence.py
```

The script preserves Open Targets query provenance, uses the existing frozen
81-gene accessible intersection as input, and does not alter the upstream
three-source convergence result.
