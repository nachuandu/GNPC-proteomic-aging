# GNPC proteomic aging analyses

Code and input-format documentation for proteomic aging analyses in GNPC, ST16 and IADRC. **No study data are distributed.** Authorized users must obtain data from the respective custodians and run these scripts in an approved environment.

## Analysis coverage

| Stage | Entry point | Scope |
| --- | --- | --- |
| GNPC preprocessing | `code/preprocess_gnpc.R` | Missing-value filtering, recovered log transform, impute.knn and z-scores |
| Protein-module clustering | `code/cluster_features.py` | Historical age-bin medians, interpolation and direct k-means; no PCA clustering |
| K comparison | `code/compare_k_resolutions.py` | K=2–25 descriptive module–age comparison |
| Trajectory | `code/fit_phate.py`, `code/principal_curve_gnpc.R` | Explicit PHATE settings and recovered R principal-curve adapter |
| Deviation score | `code/age_deviation.py` | Monotonic age GAM; spline count and penalty must be supplied |
| Protein associations | `code/protein_associations.py` | Global Spearman associations and Class I/II; within-age GAM explained-deviance R² plus Pearson P/BH |
| Wave membership | `code/define_waves.py` | Assay screening, exact gene-label aggregation and all observed membership patterns |
| ST16 comparison | `code/replicate_st16.py`, `code/st16_wave_summary.py` | Explicit assay matching, directional concordance and frozen GNPC Wave summaries |
| IADRC longitudinal analysis | `code/prepare_iadrc_inputs.py`, `code/run_single_visit.py` | Fixed GNPC features/scaling; one IADRC visit per joint refit; percentile persistence and transitions |
| Optional pathway utility | `code/reactome_ora.py` | Offline ORA with an explicit GMT release and background; not verified against historical Enrichr results |

## Release status

This is a curated code deposit, **not yet a verified end-to-end reproduction of every manuscript figure**. It combines recovered numerical routines with newly written portable adapters. See [PROVENANCE.md](PROVENANCE.md) for the distinction. Exact upstream cohort extraction, the final historical GNPC/ST16 age-GAM settings, cross-cohort annotation mapping and final pathway/figure provenance need to be checked against the saved analysis before a publication release. Do not substitute a newly chosen parameter set for the saved manuscript results.

The current IADRC entry point is `run_single_visit.py`: each query fit contains the same 1,232 GNPC participants and **one** IADRC visit. The older `run_incremental.py` and `within_age_sensitivity.py` remain for historical traceability and are superseded for the final longitudinal analysis.

## Running the code

Use an existing compatible environment. Python dependencies are in `requirements.txt`; the two R adapters require Bioconductor **impute** and CRAN **princurve**, respectively. Historical versions and current checks are distinguished in [RELEASE_CHECKS.md](RELEASE_CHECKS.md).

```bash
python -m unittest discover -s code -p 'test_*.py' -v
```

These tests use generated synthetic arrays, not cohort records. See [code/README.md](code/README.md) for stage-by-stage commands and [INPUT_FORMATS.md](INPUT_FORMATS.md) for schemas. Outputs, caches, fitted objects and all individual scores remain private. Do not upload generated outputs simply because a filename does not begin with `PRIVATE_`.

## Interpretation

- The proteomic aging score is a relative trajectory position, not age in years or an estimated aging rate.
- IADRC uses independent single-visit **joint refits**, not frozen-model prediction. Longitudinal tests use age-local GNPC-reference percentiles. Raw 0–1 scores from different refits do not constitute a common longitudinal measurement scale.
- GAM explained-deviance R² and Pearson P values answer different questions; the code keeps them separate.
- Waves are four selected continuous membership patterns. Other passing patterns are retained in the output. Threshold crossing alone does not establish biological onset or sequential recruitment.
- K-resolution comparison is descriptive, not a clustering-stability test.

No license or archived DOI is assigned by this working deposit. Author/contributor attribution and the intended license should be finalized before public release. Data access rights are separate from code access.
