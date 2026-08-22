# MBzP–CRC Phase 2A NHANES sample audit

Run timestamp (UTC): `2026-08-22T07:09:11.034213+00:00`  
Scope: sample/event audit only; no regression was performed.

## Primary audit results

- MBzP + CRC outcome available: **16,768** participants.
- Combined CRC cases: **115** unweighted cases.
- Colon cases: 111; rectal cases: 7; both type codes: 3.

The combined CRC count is the union of MCQ230A-D code 16 (Colon) and code 31 (Rectum), so colon and rectal subtype counts are not added to obtain the primary event count.

## Cycle-level CRC cases

| Cycle | MBzP + outcome N | CRC cases |
|---|---:|---:|
| 1999-2000 | 1,461 | 9 |
| 2001-2002 | 1,647 | 15 |
| 2003-2004 | 1,533 | 15 |
| 2005-2006 | 1,490 | 6 |
| 2007-2008 | 1,814 | 15 |
| 2009-2010 | 1,914 | 16 |
| 2011-2012 | 1,705 | 11 |
| 2013-2014 | 1,814 | 7 |
| 2015-2016 | 1,690 | 11 |
| 2017-2018 | 1,700 | 10 |

All values are unweighted counts. NHANES files and sampling weights are retained locally and are not committed; the manifest records the exact file map, hashes, variables, and weight columns. Phase 2B models remain deferred pending acceptance of this audit.
