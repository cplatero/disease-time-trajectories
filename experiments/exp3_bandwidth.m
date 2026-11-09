function exp3_bandwidth(mode)
%EXP3_BANDWIDTH  Permutation tests of the kinematics and bandwidth choice.
%
%   exp3_bandwidth('real')
%
%   For every cohort, it runs bw_sensitivity with the local-linear
%   sw-LOESS on the bandwidth grid (2, 3, 4, 6, 8 years) and with the
%   local-constant estimator at bw = 3. Comparisons: A+ vs A- (MEM, EF,
%   WMH) and V+ vs V- (MEM, EF, Centiloids), in both transitions.
%
%   Paper items (check numbers against the final PDF):
%     p-values of Table 4 (column p)          bw_sensitivity_<C>_ll.csv, bw = 3
%     bandwidth sensitivity (Tables S2-S3)       bw_pvalues_<C>_ll.tex
%     cross-validation (Tables S4-S5)            bw_cv_<C>_ll.tex
%     local-constant estimator (Table S9, Const.) bw_sensitivity_<C>.csv
%   It also writes calib_<C>.csv, the design of every comparison, which
%   the calibrated simulation (exp4) reproduces.
%
%   Outputs: results/<mode>/exp3_bandwidth/
%
%   Authorship: Carlos Platero (UPM), with the collaboration of Claude
%   (Anthropic).

if nargin < 1, mode = 'real'; end
C = config(mode);
out = fullfile(C.resultsDir, 'exp3_bandwidth');
common = {'nPerm', C.N_PERM, 'seed', C.SEED, 'outDir', out, 'analyses', C.ANALYSES, ...
          'dataDir', C.dataDir, 'tauDir', C.tauDir, 'wmhThreshold', C.WMH_THRESHOLD};
for c = 1:numel(C.cohorts)
    bw_sensitivity(C.cohorts{c}, 'estimator', 'linear', 'bwGrid', C.BW_GRID, ...
                   'exportCalib', true, common{:});
    bw_sensitivity(C.cohorts{c}, 'estimator', 'const', 'bwGrid', C.BW, common{:});
end
end
