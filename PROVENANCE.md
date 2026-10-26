# Code provenance and limits

## Curated update: 26 October 2026

The update adds portable code for the main computational stages. Reconstructed adapters are clearly separated from numerical routines recovered from completed analyses. No participant-level files, selected assay lists, fitted parameters or source chat transcripts are distributed.

| Component | Evidence and release status |
| --- | --- |
| Historical clustering | Recovered age-bin/interpolation/direct-k-means workflow matched the original numbered clusters and all 227 selected assays in the earlier audit. Portable command wrappers still require an authorized-input rerun. No PCA or additional curve standardization is inserted. |
| Preprocessing | Portable R adapter of a recovered historical script: -1 to NA; >50% missing exclusion; log2(RFU+1); transposed impute.knn(k=5); sample-SD scaling. Full matrix equivalence and cohort/QC extraction are not established by this deposit. |
| GNPC trajectory | PHATE adapter plus an R princurve adapter recovered from historical code. The historical early orientation was negative with age; the adapter exports the positive orientation used downstream. This source is not proof that every historical figure used that exact backend. Recheck against the saved GNPC fit before replacing archived scores. |
| Age deviations | Newly written monotonic-GAM interface. Spline count and penalty are required because historical GNPC/ST16 chats contain different trial settings; no trial is silently declared final. |
| Within-age associations | Newly written direct pyGAM implementation of the verified GAM explained-deviance/Pearson/BH procedure. The separate earlier audit established that OLS R² was not equivalent to the saved GAM results. |
| Wave aggregation | Newly written adapter preserving assay-first screening and exact-label OR aggregation; retains all patterns rather than treating the four waves as exhaustive. |
| ST16 | Newly written concordance and frozen-wave summaries. Explicit one-to-one assay mapping is required. Verify the supplied mapping and cohort-specific score provenance before interpreting count differences. |
| IADRC final | `run_single_visit.py` adapted from the completed 23 October single-visit runner. Helper paths now resolve to the code directory; the shared backend is extracted from the working fit routine. Helper hashes are included in cache fingerprints. Numerical fit and percentile operations are retained. |
| IADRC preparation | Newly written explicit-ID adapter; requires the approved eligibility/history tables and frozen GNPC scaler. It uses the completed IADRC run's log2(RFU), without +1. This difference from the recovered GNPC preprocessing recipe is documented, not silently harmonized. |
| Pathway utility | Newly written offline hypergeometric/BH implementation requiring GMT release and universe. Not equivalent by assertion to the historical Enrichr analysis. Final enrichment and network display provenance remain to be reconciled. |

## Historical clustering details retained

The recovered protocol uses one-year bins from 20 to 30, two-year bins starting at 30 through 44, and one-year bins starting at 45 through 89. Thus 44–46 overlaps 45–46, and age exactly 90 passes the eligibility range but is outside every half-open bin. Bins require at least two samples. Missing medians use NumPy interpolation with endpoint extension; fewer than two valid bins use the historical zero-fill fallback. These details are preserved for reproduction, not recommended as general binning rules.

K-means uses seed 0, n_init=10 and Lloyd's algorithm. Numbered labels depend on exact inputs, their order and software. Module selection/comparison uses the strongest positive and negative module–age correlations, not a claim of clustering stability.

## Separate curve and age-reference implementations

The R `princurve` adapter and the Python custom six-knot principal curve are separate implementations. The Python helper preserves native coordinate distances and uses PCA only to initialize the one-dimensional curve; this is unrelated to PCA before protein k-means, which is not used. Neither implementation is the full branching Slingshot algorithm.

The IADRC age reference is an **unconstrained** natural cubic smoothing spline with total EDF 6, not the monotonic pyGAM used for the primary age-deviation calculation. Each single-visit refit has its own GNPC reference curve and distribution. Raw scores from different refits are not pooled for longitudinal inference.

## Historical alternatives

The original repository deposited cumulative participant-addition code. `run_incremental.py` and `within_age_sensitivity.py` are retained as superseded historical analyses. They are not the final single-visit protocol and should not be used to describe the current Fig. 5/Supplementary sensitivity analysis.

## Validation and remaining release work

See RELEASE_CHECKS.md for checks performed on this deposit. Synthetic tests do not establish full-cohort numerical reproduction. Complete upstream cohort extraction, exact final GNPC/ST16 configurations, the historical matching table, enrichment background/release and final figure-generation scripts still need provenance reconciliation. This deposit does not claim to reproduce all assembled PowerPoint panels or every supplementary figure.

No new open-source license is assigned. The authors should confirm contributor attribution, ownership of contributed code and intended licensing before publication. Code access does not grant access to restricted data.
