# Running the analysis stages

All paths below are examples pointing to files in an approved private data environment. They are not bundled. Use **new output directories**; most adapters refuse to overwrite existing results. The code never automatically uploads results or calls data services.

## 1. Preprocess the approved GNPC control matrix

```bash
Rscript code/preprocess_gnpc.R /approved/raw_controls.rds /approved/preprocessed
```

The input must already have the approved subjects, plasma samples, QC decisions and assay columns. The recovered recipe converts -1 to NA, excludes assays with >50% missingness, uses log2(RFU+1), imputes the transposed matrix with `impute.knn(k=5)` and computes per-assay sample-SD z-scores. This adapter has not been rerun against the archived complete input matrix. Do not use it to re-log already transformed inputs.

## 2. Cluster proteins and retain trajectory features

```bash
python code/cluster_features.py --matrix /approved/standardized.csv --metadata /approved/metadata.csv --k 20 --seed 0 --clusters 9,19 --output /approved/clustering
python code/compare_k_resolutions.py --matrix /approved/standardized.csv --ages /approved/metadata.csv --output /approved/k_comparison
```

The numbered clusters are meaningful only for the recovered original ordering, preprocessing, seed and software. The historical target is 78+149=227 features. Never force a different run to this count by relabeling or truncating. Omit `--clusters` to select the most positively and negatively age-correlated modules in an exploratory rerun; that does not automatically reproduce the historical feature list.

## 3. Derive trajectory scores and age deviations

```bash
python code/fit_phate.py --matrix /approved/standardized.csv --metadata /approved/metadata.csv --features /approved/selected_features.csv --seed "${PHATE_SEED:?Set from the corresponding saved run}" --output /approved/embedding
Rscript code/principal_curve_gnpc.R /approved/embedding/PRIVATE_embedding.csv /approved/curve
python code/age_deviation.py --scores /approved/curve/PRIVATE_scores.csv --n-splines "${AGE_N_SPLINES:?Set from the final cohort-specific run}" --lam "${AGE_LAM:?Set from the final cohort-specific run}" --output /approved/deviations
```

The R curve adapter is recovered from historical GNPC code. The later IADRC implementation uses the separately provided Python six-knot curve. They are not interchangeable. For exact downstream reproduction, use the saved, positively oriented scores and residuals until the upstream run has been checked. The deviation is the raw residual, not a residual z-score.

ST16 is scored separately using its cohort-specific standardized matrix and the approved shared-feature list. Do not pool GNPC and ST16 in this step. Determine final ST16 age-GAM settings from the saved result provenance; historical chats contain multiple trials.

## 4. Protein associations and Wave membership

```bash
python code/protein_associations.py global --matrix /approved/standardized.csv --metadata /approved/scores.csv --output /approved/global_associations
python code/protein_associations.py within-age --matrix /approved/standardized.csv --metadata /approved/scores.csv --output /approved/within_age
python code/define_waves.py --associations /approved/within_age/PRIVATE_associations.csv --annotation /approved/annotation.csv --output /approved/waves
```

Global correlations are Spearman correlations. Each P-value family is BH-adjusted separately across the tested assays. `delta_rho = rho_score - rho_age`; Class I uses |rho_age| <=0.2 and |delta_rho| >0.1; Class II uses |rho_age| >=0.3 and |delta_rho| >0.1. The 0.2–0.3 interval remains unclassified. These definitions do not imply that every selected association is stronger in absolute magnitude.

Within each half-open decade [50,60), [60,70), [70,80), [80,90), the historical assay-level model is `LinearGAM(s(0,n_splines=5,spline_order=3),lam=0.6)`. Explained-deviance R² comes from that GAM; two-sided P comes from Pearson correlation, followed by BH over assays in that bin. A failed fit is not silently replaced by OLS. Screening uses R² >=0.20 and q<0.05. Gene-label membership is the logical OR over passing assays with the same exact supplied label. Empty labels are excluded; compound labels are not silently split or harmonized.

Wave patterns in age order are 1111, 0111, 0011 and 0001. The output also retains 0000 and every other selected pattern. Use `--r2` for a sensitivity threshold; do not edit counts manually.

## 5. ST16 replication

Run the global and within-age commands separately on ST16, then:

```bash
python code/replicate_st16.py --gnpc /approved/GNPC_global.csv --st16 /approved/ST16_global.csv --mapping /approved/approved_assay_mapping.csv --output /approved/ST16_concordance
python code/st16_wave_summary.py --associations /approved/ST16_within_age.csv --annotation /approved/ST16_annotation.csv --gnpc-waves /approved/waves/PRIVATE_gene_label_membership.csv --output /approved/ST16_waves
```

The mapping adapter retains only one-to-one pairs and excludes ambiguity explicitly. Directional concordance is sign agreement of delta_rho, not independent significance replication. Zero differences are reported separately. The denominator is evaluable matched pairs, not all GNPC assays. Exact historical mapping equivalence still needs checking with the approved local mapping. ST16 Wave summaries retain the GNPC definitions; they do not select new ST16 waves.

## 6. IADRC single-visit analysis

Prepare the input bundle using `prepare_iadrc_inputs.py --help`, or use the already checked private bundle described in INPUT_FORMATS.md.

```bash
python code/run_single_visit.py --input /approved/iadrc_bundle --output /approved/iadrc_single_visit --workers 4 --threads 2
```

A `--limit 1` smoke run performs the GNPC reference fit and one query fit; it is not the final analysis. The completed run uses 443 query fits plus one GNPC-only fit, and summarizes 133 baseline-to-last pairs. Each query has exactly 1,233 rows. The script checks the 134-person all-annotated-visits-healthy inclusion rule, fixed inputs, unique IDs, convergence and cache fingerprints. Checkpoints allow resuming the same protocol. Age percentiles use the within-fit GNPC age residual distribution, Gaussian age weights (five-year bandwidth) and a weighted mid-distribution. Cutoffs are 100/3 and 200/3. The GNPC sensitivity diagnostic excludes each queried GNPC person from their own reference distribution.

`run_incremental.py` and `within_age_sensitivity.py` are historical alternatives, not this final method.

## 7. Optional offline pathway analysis

```bash
python code/reactome_ora.py --membership /approved/waves/PRIVATE_gene_label_membership.csv --background /approved/tested_background.csv --gmt /approved/reactome.gmt --release "${REACTOME_RELEASE:?Supply the exact GMT release}" --output /approved/reactome_rerun
```

This is a newly written explicit-background utility, not a verified recreation of the manuscript's historical Enrichr calculations. It tests all pathways with background overlap and applies BH within each age bin; it does not reproduce the final displayed-pathway selection or network layout.
