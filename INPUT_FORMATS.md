# Input formats, without data

This document describes interfaces only. No actual rows, identifiers, numerical arrays, assay lists, input-file checksums or study outputs are included. Fields named `pseudotime` in the code mean the trajectory-derived proteomic aging score. Fields named `gap` mean the proteomic deviation score.

## General rules

- CSV text is UTF-8 with a header row; NumPy input arrays are numeric `.npy` files loaded with `allow_pickle=False`.
- A sample is one protein measurement visit; a participant may have multiple samples.
- Ages are numeric years at protein sampling. Keep within-person visit order explicit.
- Do not infer negative clinical findings from missing or unknown flags.
- Array row order must match the corresponding metadata CSV exactly. Feature columns must have the same assay identity and order in GNPC and IADRC; dimensions alone are not evidence of matching.
- Standardization parameters and the selected assay order are prepared locally from GNPC. They are not estimated anew from IADRC by these scripts.

## Cumulative-addition inputs

The runner expects the following files in a user-supplied approved input directory. The strict dimensions below describe this particular archived analysis, not a general-purpose interface.

| File | Required form |
|---|---|
| `GNPC_X.npy` | Finite numeric array, 1,232 GNPC rows by 227 fixed assay features, in GNPC reference-standardized units. |
| `IADRC_reference_scaled.npy` | Finite numeric array, 651 archived candidate visit rows by the same 227 features, standardized using the fixed GNPC parameters. |
| `GNPC_reference.csv` | One row per GNPC array row. Required columns: `sample_id` (unique key), `age` (numeric years), `pseudotime` (numeric saved historical reference score used only for comparison). |
| `IADRC_visit_metadata.csv` | One row per archived IADRC array row. Required columns: `sample_id` (unique key), `participant` (participant key), `visit` (integer order), `age` (numeric years), `flag_interpretation` (clinical eligibility category). |
| `people.csv` | One row per candidate participant. Required columns: `participant`, `original_240`, `followup_status`, `all_record_status`. These are eligibility annotations, not inferred diagnoses. |
| `visits.csv` | Clinical annotation history. Required columns: `participant`, `healthy_current` (Boolean or a Boolean string). This file must include the earlier annotations used to exclude people with earlier abnormal flags. |
| `input_manifest.json` | JSON object with an `input_sha256` mapping from local input basenames to SHA-256 values. Generate this locally; never commit the real manifest. |

Eligibility labels expected by the archived runner are `original_240 == "Yes"`, `followup_status == "Later all healthy flags"`, and `all_record_status == "All recorded flags healthy"`. The selected visits must have `flag_interpretation == "H"`. The runner asserts 134 selected participants and 443 selected protein visits. Selected baseline visit order must equal 1. It expects 133 baseline-to-last pairs with follow-up.

The 651 candidate rows are read for selection and alignment only. Each refit includes GNPC plus visits of the cumulative subset of the selected 134 participants; visits from the other candidate participants do not enter a refit. Missing annotation is not a negative flag. Correct upstream clinical adjudication and historical linkage are required; these scripts do not reconstruct them from raw clinical records.

`participant` values must have the same type across the eligibility and visit tables. The exact historical participant keys and their sorting order matter for reproducing the prespecified random addition order. Changing or reassigning keys may change intermediate refits even when the random seed is unchanged.

## Input-integrity manifest

The manifest can be described structurally as an object containing the key `input_sha256`, whose value is an object of `local_basename: lowercase_sha256_string` entries. Include all required data inputs. Do not use this structure to publish actual file hashes or paths.

Numerical helper modules are loaded from `code/` in this portable deposit rather than copied into the input directory. They therefore do not need to be listed in the data manifest.

## Inputs for the K-resolution comparison

| Input | Required form |
|---|---|
| Expression CSV | Samples by assay features. First column is the unique sample key used as the row index; subsequent columns are numeric standardized protein measurements. No missing values at this stage. |
| Age CSV | One row per sample, with columns `sample_id` and `age`; the keys must align with the expression index. |

This comparison begins with the already processed GNPC matrix, not raw RFU. Upstream missing-value handling, log transformation and standardization must be reproduced separately. The script does not invent or repeat those steps.

## Locally generated outputs

The analysis writes scores, models, participant orders, diagnostics and result tables into the output path supplied by the user. These are **not repository contents** and must remain local or in the approved analysis environment. For the additional within-age sensitivity script, provide the cumulative output directory and its `PRIVATE_GNPC_final.csv` metadata output.

## Not yet deposited

Raw-provider file schemas, the complete preprocessing pipeline, ST16 replication inputs, Reactome input schemas and exact figure-source mapping still need a separate code-and-documentation review. This initial document should not be interpreted as covering those interfaces.
