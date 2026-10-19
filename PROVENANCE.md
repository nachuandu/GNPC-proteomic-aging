# Code provenance and release boundaries

## Initial deposit

The cumulative IADRC runner is adapted from the completed working cumulative-addition script. Its numerical operations are unchanged. The three helper-module paths now resolve relative to the deposited code directory instead of the restricted input directory.

The numerical helper files contain only the functions used by that runner, extracted from the working implementations. Historical plotting routines, machine-specific paths, unrelated exploratory analyses and study-data checksums have been omitted. These helpers implement a custom single principal curve and an unconstrained fixed-effective-degrees-of-freedom natural cubic smoothing spline; they are not the full Slingshot package or a monotonic GAM.

The separate within-age sensitivity script extracts the existing leave-one-reference-person-out weighted-percentile diagnostic and accepts explicit input/output paths. It compares each cumulative fit with the GNPC-only fit.

The K-resolution script is a portable transcription of the recovered baseline age-bin and k-means numerical steps, followed by the same increasing/decreasing module-age correlation comparison. It is not a stability test. Its dataset-wide equivalence still requires an authorized-input rerun before a publication release.

## Historical binning details retained for reproducibility

The recovered protocol uses 1-year bins from 20 to 30, 2-year bins starting at 30 through 44, and 1-year bins starting at 45 through 89. Consequently, 44-46 overlaps 45-46, and an age of exactly 90 passes the eligibility range but is outside every half-open bin. These are retained historical implementation details, not recommended general bin definitions.

Bin medians require at least two people. Missing bins use NumPy interpolation, including its endpoint extension behavior; curves with fewer than two valid bins follow the historical zero-fill fallback. No PCA, extra per-curve z-score or pre-clustering smoothing is applied by this comparison. Display smoothing is not part of the clustering input.

## Testing and limits

The included tests use only in-memory synthetic arrays and check numerical identities and invariances. Five checks passed in an existing local environment using Python 3.13.5, NumPy 2.4.2, SciPy 1.17.0, scikit-learn 1.9.0 and pandas 3.0.0. This is distinct from the archived analysis environment recorded in `requirements.txt`; a new clean installation was not tested. Compilation and synthetic checks are not a full-cohort reproduction or a verification of manuscript result values. No study data, fitted parameters, selected assay list, results or exact software-environment archive are included.

This initial repository does not yet contain a complete end-to-end reproduction of the manuscript. Preprocessing, all ST16 analyses, protein associations, pathway enrichment and complete figure-generation code remain to be curated and tested.

## Attribution and licensing

No new open-source license is assigned in this initial private repository. Before publication, the authors should confirm contributor attribution, ownership of contributed principal-curve code, compatibility with dependencies and the intended software license. Access to this repository does not grant rights to the underlying cohort data.
