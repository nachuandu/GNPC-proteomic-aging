# GNPC proteomic aging analyses

Analysis code and input-format documentation for a study of proteomic aging in GNPC, ST16 and IADRC.

**Code only. No study data are included.** This includes no raw or processed measurements, participant records, assay lists, fitted models, individual scores, result tables, manuscript files or figure images.

## Status

This is an initial private working repository, not a complete publication release. The first code deposit covers the IADRC cumulative participant-addition analysis, its numerical dependencies and a retrospective K-resolution comparison. Do not yet describe this repository as a complete reproduction of every main and supplementary figure.

The repository does not develop or distribute a replacement for PHATE or the full Slingshot algorithm. The supplied custom single principal curve has no lineage tree or branching step.

## Contents

- `code/`: analysis scripts and numerical helper functions; see its README for scope and commands.
- `INPUT_FORMATS.md`: required input filenames, dimensions, fields, types and alignment rules. No real records or example data files.
- `DATA_POLICY.md`: data exclusions and access boundaries.
- `PROVENANCE.md`: the relationship between the deposited code and the working analysis, including limitations.
- `requirements.txt`: versions recorded in the completed cumulative-addition run, plus a dependency whose version still needs checking.
- `.gitignore`: blocks data, outputs, models, logs, local environments and common secret files. It is not a substitute for manual review.

## Getting started

Use an existing compatible Python environment. The recorded analysis used Python 3.12.6. Installation of dependencies is the user's choice; no environment is bundled here.

Run the code-only synthetic numerical checks without providing any cohort data:

```bash
python code/test_core.py
```

These tests generate arrays in memory. They do not distribute cohort data or establish reproduction of the manuscript's numerical results.

Authorized users must obtain cohort inputs through the relevant data custodians, prepare them in their approved environment, and supply their own paths. See `INPUT_FORMATS.md` and `code/README.md`.

## Interpretation

- Cumulative addition tests sensitivity of the GNPC reference axis to IADRC inclusion. The final analysis is a joint fit, not independent prediction from a frozen GNPC model.
- Within-age percentiles are relative to age-local GNPC residuals. Baseline tertiles are kept fixed during follow-up.
- Module-age correlations used to compare K are descriptive, post-selection summaries. They do not establish clustering stability or independent predictive performance.
- A proteomic aging score is a relative trajectory position, not chronological years or an independently calibrated biological age.

## Before a publication release

Complete the remaining preprocessing, ST16, protein-association, pathway and figure-code inventory; check the exact environment and figure mapping; confirm code attribution and licensing with the authors; and test the complete pipeline in the approved data environment. No open-source license or DOI is assigned by this initial deposit.
