# Private input schemas

Inputs and outputs remain in the approved data environment. None of the files described below are distributed. Column names describe data structure only.

## Shared table conventions

- CSV uses UTF-8 and headers. `sample_id` is a string; preserve leading zeros.
- Processed expression matrices: first column `sample_id`, then unique assay-feature columns, one row per sample. All processed values must be finite.
- Metadata: unique `sample_id`, numeric `age`, and, where needed, `proteomic_aging_score` and `proteomic_deviation_score`. Align by explicit sample ID, never by unverified row position.
- Feature selection: CSV with one unique `feature` per row, in the intended model input order. The actual selected list is not bundled.
- Annotations: unique `feature`, exact `gene_label`. Use empty labels explicitly for unannotated assays. Supply the original approved labels; the code does not infer aliases or split compound targets.
- Global associations: `feature`, `rho_age`, `rho_score`, `delta_rho`, `p_age`, `q_age`, `p_score`, `q_score`, `class`.
- Within-age associations: one row per feature per bin, `feature`, `age_bin`, `R2_deviation`, `pearson_r`, `p_pearson`, `q_pearson`, `R2_age`. Age bins are `50_60`, `60_70`, `70_80`, `80_90`.
- ST16 mapping: `gnpc_feature`, `st16_feature`; ambiguous and many-to-one pairs are not silently accepted.

## GNPC preprocessing

The RDS input is a numeric matrix with sample row names and assay column names, already restricted to approved plasma controls and QC-passing samples. Raw missing sentinel -1 is converted to NA. Do not pass an already log-transformed matrix to the RFU adapter. The upstream clinical harmonization and selection decisions are not inferred by this script.

The private scaler CSV has `feature`, `mean_log2`, `sd_log2`; SD is the sample SD used in the saved GNPC standardization. Do not re-estimate it on IADRC.

## IADRC bundle for the final single-visit runner

| File | Schema / constraint |
| --- | --- |
| `GNPC_X.npy` | Finite 1232 x 227 frozen standardized GNPC matrix |
| `IADRC_reference_scaled.npy` | Finite 651 x 227 eligible-visit IADRC matrix in the same feature order, transformed using fixed GNPC parameters |
| `GNPC_reference.csv` | Unique string `sample_id`, numeric `age`; rows aligned to GNPC_X |
| `IADRC_visit_metadata.csv` | Unique string `sample_id`, `participant`, numeric `visit`, numeric `age`, `flag_interpretation`; rows aligned to IADRC matrix |
| `people.csv` | Unique `participant`, `original_240`, `followup_status`, `all_record_status` |
| `visits.csv` | Complete reviewed clinical history: `participant`, `healthy_current` |
| `input_manifest.json` | Object `input_sha256`: filename-to-SHA256 mapping of the six inputs above |

The retained participants have `original_240="Yes"`, `followup_status="Later all healthy flags"`, and `all_record_status="All recorded flags healthy"`. All their available annotated history rows must have `healthy_current=true`; eligible proteomic visits must have `flag_interpretation="H"`. Missing or unknown flags must not be recoded as healthy. The supplied clinical curation is responsible for completeness of history, chronology and flag semantics. Participant identifiers must use the same type and representation in all three IADRC tables.

The runner checks 134 people, 443 selected visits and 133 baseline-to-last pairs. Within each person, `visit` is the original eligible-visit ordering; selected baseline must be 1. Use the actual original ordering, not an order inferred from file rows. GNPC and IADRC sample IDs must be disjoint.

`prepare_iadrc_inputs.py` can assemble this bundle from approved processed inputs. Supply its required paths shown by `--help`: standardized GNPC matrix and metadata, positive-RFU IADRC matrix, frozen feature list/scaler, visit metadata, participant selection table and complete reviewed history. It uses explicit sample-ID joins and then invokes the final runner's inclusion checks. It does not recover identities or diagnoses from expression values.

## Optional pathway inputs

- Background CSV: `gene_label`, the actual tested/eligible universe for the specified analysis.
- GMT: tab-separated pathway ID/name, description, then gene labels. Provide its exact release identifier.
- Membership: the generated gene-label table with four binary age-bin columns.

No gene lists are sent to a remote service. The offline utility is not yet verified as a reconstruction of the manuscript's historical Enrichr results.
