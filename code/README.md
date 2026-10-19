# Analysis code

This initial code deposit is incomplete at manuscript level. All commands below require a compatible Python environment. Cohort inputs and all outputs must remain outside the repository in an approved analysis environment.

## Included scripts

| Script | Scope |
|---|---|
| `run_incremental.py` | Participant-by-participant cumulative GNPC-IADRC refitting, final persistence analysis and global reference-score diagnostics. |
| `principal_curve_methods.py` | Native-coordinate single principal curve with six natural-spline knots by default. |
| `GAM_methods.py` | Integrated squared-second-derivative penalty matrix. |
| `fixed_GAM_methods.py` | Unconstrained GNPC-only natural cubic smoothing spline with fixed total EDF. |
| `within_age_sensitivity.py` | Additional leave-one-GNPC-person-out within-age percentile sensitivity across cumulative fits. |
| `compare_k_resolutions.py` | Retrospective K=2-25 comparison of increasing and decreasing module-age correlations under recovered binning rules. |
| `test_core.py` | Code-only synthetic numerical checks; no cohort inputs or stored example data. |

## Cumulative joint trajectory

```bash
python code/run_incremental.py --input /approved/input --output /approved/output --threads 8
```

The illustrative `/approved/...` locations are placeholders, not actual paths or access grants. The default final step is 134, giving 135 fits including GNPC-only. The random participant order is fixed with seed 20261018; PHATE uses seed 20261012, `knn=50`, `t=50`, metric SMACOF, and the settings embedded in the runner. Each participant's visits enter together. All longitudinal results come from the final fit, not mixed addition-step scores.

`--seed-controls` additionally runs the archived endpoint seed checks. `--last-step` can stop the audit early, but such a partial run does not generate the full final longitudinal result.

The output directory must be fresh for a new protocol. The archived runner resumes saved steps; reusing a directory from another seed, dataset, helper version or parameter configuration is unsafe. The input manifest verifies inputs, but does not make mismatched checkpoints valid.

## Age-local percentile sensitivity

After the complete cumulative run:

```bash
python code/within_age_sensitivity.py --results /approved/output --output /approved/diagnostics
```

The script reads 135 saved steps and final GNPC metadata, refits the GNPC-only age-score spline for each step, and excludes each queried GNPC individual from that step's weighted reference CDF. The mean and 95th percentile of absolute shifts are descriptive summaries, not confidence intervals.

## K-resolution comparison

```bash
python code/compare_k_resolutions.py --matrix /approved/expression.csv --ages /approved/ages.csv --output /approved/k_comparison
```

The matrix is already processed and standardized; this script does not start from raw RFU. It uses the historical piecewise bins, minimum bin size 2, interpolation, no PCA and no pre-clustering smoothing, then k-means with seed 0 and `n_init=10`. It selects the strongest positive and strongest negative individual-level mean-module-age Spearman correlations at each K, and reports their mean absolute value. A missing correlation direction is reported rather than replaced by an arbitrary module. Cluster labels have no intrinsic biological identity across different K values.

See `PROVENANCE.md` for historical overlapping-bin and endpoint details. This score is a descriptive age-association criterion, not clustering stability, cross-validation or independent evidence that a given K is optimal.

## Security and scope

No script transmits data or downloads restricted data. Outputs may contain sensitive records and must not be uploaded. Input formats are described in `../INPUT_FORMATS.md`. Exact manuscript-result reproduction and all missing analysis stages still require a separate authorized-environment audit.
