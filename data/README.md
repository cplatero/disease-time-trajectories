# Input data

The ADNI and NACC data are not distributed with this repository: their Data Use Agreements do not allow it. This folder holds the input files on the machine where the experiments are run; `*.csv` files here are ignored by git.

`run_all('demo')` creates synthetic files with the same layout in `data/demo/`.

## Visit tables: `<COHORT>_<TRANSITION>.csv`

`COHORT` is `ADNI` or `NACC`; `TRANSITION` is `CU2MCI` (cognitively unimpaired to MCI) or `MCI2Dem` (MCI to dementia). One row per visit.

| Column | Meaning |
|---|---|
| `ID` | subject identifier (ADNI RID, NACC ID) |
| `AGE_visit` | age at the visit, years |
| `ClinGroup` | 1 = stable, 2 = converter during follow-up, 3 = already in the next stage at baseline; empty if the subject is not used in the clinical row (e.g. ADNI SMC) |
| `MEM`, `EF` | harmonized memory and executive function composites (ADNI: ADNI_MEM, ADNI_EF; NACC: ADSP-PHC PHC_MEM, PHC_EXF) |
| `CENTILOIDS` | amyloid PET, Centiloids |
| `WMH` | white matter hyperintensity volume, ml |
| `A_bl` | baseline amyloid status, 1 = A+, 0 = A−, empty without PET (ADNI: PET status; NACC: `CENTILOIDS_bl` > 20) |
| `CENTILOIDS_bl`, `WMH_bl` | first valid value of each subject |
| `onset_age_clin` | converters only: age at the first visit with the new diagnosis |

Subjects of the CU→MCI file: cognitively unimpaired (stable or converter) and stable MCI. Subjects of the MCI→dementia file: MCI (stable or converter) and dementia.

## Transition-age estimates: `tau/tau_<COHORT>_<TRANSITION>.csv`

One row per subject: `ID`, `tau_1`, …, `tau_10`, the transition ages (years) estimated by the 10 Leaspy models of the disease progression model. Empty or `Inf` values are ignored. `load_cohort` uses their mean.

## Building the files

The author built them with `prepare/run_prepare.m`, which calls

- `export_ADNI_csv.m`: ADNIMERGE, UW neuropsychological summary scores, amyloid PET (Centiloids) and UC Davis WMH tables, through the selection of `prepare/aux_DPM/prepare_ADNI_Amyloid.m`;
- `export_NACC_csv.m`: the NACC longitudinal tables of each transition (UDS with ADSP-PHC harmonized cognition, amyloid PET and WMH);
- `export_tau_csv.m`: the transition-age estimates of the DPMs.
