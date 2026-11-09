# Comparing subgroup trajectories on a model-estimated disease time: subject-weighted local-linear smoothing and simulation-calibrated permutation tests

MATLAB code of the paper

> C. Platero. *Comparing subgroup trajectories on a model-estimated disease time: subject-weighted local-linear smoothing and simulation-calibrated permutation tests.* Submitted to *Computer Methods and Programs in Biomedicine*.

Longitudinal studies often compare a marker between subgroups (for example, amyloid-positive vs. amyloid-negative subjects) along age or follow-up time. A disease time estimated with a disease progression model (DPM) reflects the stage of the disease better, but comparing subgroups on this axis raises four problems: unequal numbers of visits per subject, subgroups that occupy different regions of the axis, dependent visits within each subject, and error in the estimated transition age of each subject.

The framework has four components:

1. **Disease-time axis.** Each visit is placed at `t* = age at visit − τ̂ᵢ`, where `τ̂ᵢ` is the transition age estimated by the DPM (mean of 10 Leaspy models).
2. **sw-LOESS.** A non-parametric, subject-weighted local-linear smoother. Each subject contributes one local weight, whatever his number of visits, and no global shape is imposed on the curve.
3. **Permutation test.** Subject labels are permuted and the mean absolute difference between the two curves, `T`, is computed on their common range.
4. **Calibrated p-value.** Null replicates reproduce the design of each comparison and the transition-age error of its cohort: `p_cal = (1 + #{T_null ≥ T_obs}) / (1 + R)`.

## Quick start (no data needed)

```matlab
cd disease-time-trajectories
setup
run_all('demo')                 % every experiment on synthetic data (about 10 minutes)
edit tutorial/tutorial_framework.m   % step-by-step walk-through of one comparison
```

The synthetic cohorts of `data/demo/` have the same layout as the real ones, but their numbers have no clinical meaning.

## Requirements

- MATLAB R2020a or later, with the Statistics and Machine Learning Toolbox (`prctile`, `iqr`).
- `sim_framework.m` (experiment 1) also runs in GNU Octave 8.

## Paper items and experiments

Each experiment is a function in `experiments/` that writes its outputs to `results/<mode>/<experiment>/`. `run_all(mode)` runs them in order.

| # | Experiment | Paper items | Main outputs | Data |
|---|---|---|---|---|
| 1 | `exp1_generic_simulation` | Fig. 2, Table 3; Tables S7–S8 (type I error); bias and coverage tables (not in the paper) | `typeI/sim_typeI.png`, `typeI/sim_typeI.tex`, `power/sim_power.tex`, `power/sim_bias.tex`, `power/sim_coverage.tex` | none |
| 2 | `exp2_transition_age_error` | Table 1 | `tau_error_table.tex/.csv` | ADNI, NACC |
| 3 | `exp3_bandwidth` | column *p* of Table 4; Tables S2–S5 (bandwidth sensitivity and cross-validation) | `bw_sensitivity_<C>_ll.csv`, `bw_pvalues_<C>_ll.tex`, `bw_cv_<C>_ll.tex`, `calib_<C>.csv` | ADNI, NACC |
| 4 | `exp4_calibrated_pvalues` | column *p*<sub>cal</sub> of Table 4; Fig. 3; Table S9 | `table_p_pcal.tex/.csv`, `fig_sim_cal.png`, `sim_out_<C>/sim_calibrated.tex` | outputs of 3 |
| 5 | `exp5_kinematics` | Figs. 4–5; Tables S11–S12 | `fig/HN_<C>_<CU2MCI\|MCI2AD>_v3.png`, `metrics_table_<C>_ll_v3.tex` | ADNI, NACC + outputs of 3, 4 |
| 6 | `exp6_calibration_checks` | Table S10 (checks of the calibration) | `sim_checks_summary.tex/.csv` | outputs of 3, 4 |

`<C>` is the cohort, ADNI or NACC. Run times with the paper settings (`'real'`): experiment 1 takes several hours, 3 and 4 about one hour per cohort with `C.ANALYSES = {'A', 'V'}` and roughly twice that with the default of `config.m`, which adds the A × V comparisons of Fig. 3 and Table S9, and the rest a few minutes.

## Running on ADNI and NACC

The ADNI and NACC data cannot be redistributed. Researchers with access can rebuild the input files:

1. Request access to [ADNI](https://adni.loni.usc.edu/) and [NACC](https://naccdata.org/) (including the ADSP-PHC harmonized measures for NACC).
2. Build the visit tables described in [`data/README.md`](data/README.md) (`data/ADNI_CU2MCI.csv`, ...) and the transition-age estimates (`data/tau/tau_ADNI_CU2MCI.csv`, ...). The scripts of `prepare/` show how the author built them from the original files.
3. Run `run_all('real')`.

`config.m` holds the paths, the thresholds (A+: Centiloids > 20; V+: WMH > 10 ml), the bandwidth (3 years), the number of permutations and replicates, and the transition-age error of Table 1 used by the calibration.

## Repository structure

```
config.m            paths and settings ('real' or 'demo')
setup.m             adds the folders to the MATLAB path
run_all.m           runs the six experiments in order
experiments/        one function per group of tables and figures of the paper
src/                framework and analysis functions
  load_cohort.m       reads the CSV files and builds the disease-time axis
  sw_loess_ll.m       subject-weighted local-linear smoother
  coverage_range.m    coverage range of a subgroup
  bw_sensitivity.m    permutation tests, bandwidth analysis, design of each comparison
  plot_kinematics.m   natural-history figures and trajectory metrics
  plot_sim_cal.m      Fig. 3: type I error of each comparison and p vs. p_cal
  calib_tau_error.m   error model of the transition age
  sim_framework.m     generic and calibrated simulations
  make_demo_data.m    synthetic cohorts
tutorial/           step-by-step example of the framework
prepare/            author-side conversion of the original files to CSV
data/               input files (not versioned; see data/README.md)
results/            outputs, one folder per mode and experiment
```

## Citation

If you use this code, please cite the paper above. The disease progression models follow [Platero & Bengoa, *Comput. Methods Programs Biomed.* 274 (2025) 109162](https://doi.org/10.1016/j.cmpb.2025.109162).

## Authorship and licence

Carlos Platero (Universidad Politécnica de Madrid), with the collaboration of Claude (Anthropic). Released under the MIT licence (see `LICENSE`).
