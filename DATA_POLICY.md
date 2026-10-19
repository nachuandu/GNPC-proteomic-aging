# Data policy

## What this repository contains

Only source code, dependency information and prose descriptions of input data formats.

## What must not be committed

- GNPC, IADRC or ST16 measurements, including processed or standardized expression matrices.
- Participant-level metadata, visit histories, clinical flags, linkage keys or identifiers.
- Lists of actual participant IDs, sample IDs or selected assay IDs.
- Saved coordinates, proteomic aging scores, deviation scores, percentiles, age curves, fitted models or checkpoints.
- Result tables, reports containing calculated study results, figures, screenshots, manuscript drafts or internal data-use agreements.
- Credentials, cookies, tokens, private keys, cloud access links, personal machine paths or deployment logs.

The exclusion applies even when a file is small or described as de-identified. A private GitHub repository is not an approved storage location for the restricted study data.

## Data access

No access to data is granted by access to this code repository. Authorized users must follow the relevant provider's current application, consent, institutional approval and data-use agreement requirements.

- GNPC: consult the GNPC dataset entry and access process in the AD Discovery Portal / AD Workbench: https://www.alzheimersdata.org/discovery-portal
- IADRC: consult the Indiana Alzheimer's Disease Research Center resource-request process: https://medicine.iu.edu/research-centers/alzheimers/for-researchers/resource-request
- ST16: consult the original lifespan plasma proteomics publication and its supplementary materials, and follow the applicable reuse conditions: https://doi.org/10.1038/s41591-019-0673-2

The code does not download data or circumvent these processes. Any locally generated outputs must remain in an appropriately approved environment. The access information above is a pointer to the custodians, not a new data-sharing authorization.

## Publication readiness

This code-only repository does not by itself satisfy all manuscript Data Availability obligations. Any source-data sharing must be arranged separately with the relevant permissions; it is not authorized in this repository.
