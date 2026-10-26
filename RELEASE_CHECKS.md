# Checks on the 26 October 2026 code deposit

## Completed locally

- Parsed all 20 Python scripts successfully.
- Checked `--help` for all 10 newly supplied command-line entry points.
- Ran 13 synthetic tests: 12 passed; the pyGAM-specific test was skipped because pyGAM was not installed in the local check environment.
- Verified, by Python AST comparison, that the final IADRC runner's `selected_input`, `age_reference`, `fit_case` and `aggregate` functions match the completed working single-visit source exactly. Path resolution and helper fingerprints were adapted outside those functions.
- Reviewed the explicit upload file list: code, Markdown documentation and dependency specification only. No matrices, metadata, clinical/genetic records, feature lists, scores, fitted models, outputs or original chat transcripts included.
- Searched the deposit for machine-specific private paths and common credential markers; none found.

Checks used an existing environment: Python 3.13.5, NumPy 2.4.2, pandas 3.0.0, SciPy 1.17.0, scikit-learn 1.9.0 and threadpoolctl 3.6.0. This differs from the recorded Python 3.12.6 analysis environment; requirements.txt preserves available historical pins rather than claiming an exact newly tested lockfile.

## Not established by these checks

- The full GNPC/ST16 analysis was not rerun on authorized cohort matrices during this code-curation task.
- The 443 PHATE query refits were not repeated locally; the completed working single-visit implementation was adapted.
- R preprocessing and R princurve were not executed locally; R was unavailable in the check environment.
- No fresh clean-environment dependency installation was tested.
- Reconstructed mapping, historical age-GAM configuration, enrichment and final figure assembly are not asserted to reproduce every published-panel number.

These distinctions should remain visible when preparing the final code-availability statement and archived software release.
