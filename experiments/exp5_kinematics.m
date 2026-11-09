function exp5_kinematics(mode)
%EXP5_KINEMATICS  Natural-history figures and trajectory metrics.
%
%   exp5_kinematics('real')
%
%   Paper items (check numbers against the final PDF):
%     Figs. 4-5    fig/HN_<C>_<CU2MCI|MCI2AD>_v3.png   (rows: clinical
%                  groups, amyloid status, vascular burden; columns: MEM,
%                  EF, Centiloids, WMH)
%     Tables S11-S12 metrics_table_<C>_ll_v3.tex  (start/end values, slopes,
%                  p and p_cal of every curve)
%   The p-values come from exp3 and p_cal from exp4, so the figures and
%   Table 4 use the same tests. The p-values of the marker that defines
%   each stratification are computed here (1000 permutations).
%   Outputs: results/<mode>/exp5_kinematics/
%
%   Authorship: Carlos Platero (UPM), with the collaboration of Claude
%   (Anthropic).

if nargin < 1, mode = 'real'; end
C = config(mode);
in3 = fullfile(C.resultsDir, 'exp3_bandwidth');
in4 = fullfile(C.resultsDir, 'exp4_calibrated_pvalues');
out = fullfile(C.resultsDir, 'exp5_kinematics');
for c = 1:numel(C.cohorts)
    co = C.cohorts{c};
    pFile    = fullfile(in3, sprintf('bw_sensitivity_%s_ll.csv', co));
    pcalFile = fullfile(in4, sprintf('sim_out_%s', co), 'sim_results.csv');
    if ~exist(pFile, 'file'),    pFile = '';    warning('exp5:p', 'No exp3 output: p computed here.'); end
    if ~exist(pcalFile, 'file'), pcalFile = ''; warning('exp5:pcal', 'No exp4 output: tables without p_cal.'); end
    plot_kinematics(co, 2, C.SHOW_CI, 'pFile', pFile, 'pcalFile', pcalFile, ...
        'saveDir', fullfile(out, 'fig'), 'outDir', out, ...
        'dataDir', C.dataDir, 'tauDir', C.tauDir, 'wmhThreshold', C.WMH_THRESHOLD);
end
end
