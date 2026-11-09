function run_all(mode, steps)
%RUN_ALL  Reproduce every table and figure of the paper, in order.
%
%   run_all('demo')          % synthetic data, a few minutes, no data access needed
%   run_all('real')          % ADNI and NACC (needs the CSV files of data/)
%   run_all('real', [3 4 5]) % only some experiments
%
%   1  exp1_generic_simulation    simulation with known truth (Fig. 2, Table 3, Tables S7-S8)
%   2  exp2_transition_age_error  error of the transition age (Table 1)
%   3  exp3_bandwidth             permutation tests and bandwidth (Table 4 p, Tables S2-S5)
%   4  exp4_calibrated_pvalues    calibrated p-values (Table 4 p_cal, Fig. 3, Table S9)
%   5  exp5_kinematics            natural-history figures (Figs. 4-5, Tables S11-S12)
%   6  exp6_calibration_checks    sensitivity of the calibration (Table S10)
%
%   Experiment 1 needs no data. 3 -> 4 -> 5 and 3 -> 4 -> 6 must run in
%   this order. Results are written to results/<mode>/.
%   In 'demo' mode the synthetic data are created first (make_demo_data).
%
%   Authorship: Carlos Platero (UPM), with the collaboration of Claude
%   (Anthropic).

if nargin < 1, mode = 'demo'; end
if nargin < 2, steps = 1:6; end
C = config(mode);
if strcmp(C.mode, 'demo') && ~exist(fullfile(C.dataDir, 'ADNI_CU2MCI.csv'), 'file')
    make_demo_data(C.dataDir, C.tauDir);
end
fprintf('\n=== run_all (%s): results in %s ===\n', C.mode, C.resultsDir);
EXP = {@exp1_generic_simulation, @exp2_transition_age_error, @exp3_bandwidth, ...
       @exp4_calibrated_pvalues, @exp5_kinematics, @exp6_calibration_checks};
for k = steps(:)'
    fprintf('\n=== Experiment %d: %s ===\n', k, func2str(EXP{k}));
    t0 = tic;
    EXP{k}(mode);
    fprintf('=== Experiment %d done in %.1f min ===\n', k, toc(t0) / 60);
end
end
