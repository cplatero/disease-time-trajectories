%% Tutorial: comparing subgroup trajectories on disease time
% This script walks through the framework of the paper on ONE comparison:
% memory (MEM) in amyloid-positive (A+) versus amyloid-negative (A-)
% subjects, around the transition from cognitively unimpaired to MCI.
% It uses the synthetic cohort of data/demo, so it runs without access to
% ADNI or NACC; with config('real') it runs on the real data.
%
% Run it section by section (Ctrl+Enter in the MATLAB editor), or open it
% as a Live Script (right click > Open as Live Script).
%
% The paper's numbers are produced by the experiments of experiments/;
% this script only shows, step by step, what they compute.
%
% Authorship: Carlos Platero (UPM), with the collaboration of Claude
% (Anthropic).

%% 0. Set-up
close all;
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setup.m'));
C = config('demo');                                   % 'real' for ADNI and NACC
if strcmp(C.mode, 'demo') && ~exist(fullfile(C.dataDir, 'ADNI_CU2MCI.csv'), 'file')
    make_demo_data(C.dataDir, C.tauDir);
end
if ~exist('DO_PLOTS', 'var'), DO_PLOTS = true; end
COHORT = 'ADNI';  TRANS = 1;                          % 1 = CU -> MCI, 2 = MCI -> dementia
MARKER = 'MEM';   BW = C.BW;   T_RANGE = [-10 5];   RANGE_P = 10;   NQ = 150;

%% 1. The disease-time axis
% A disease progression model (DPM, Leaspy) gives each subject an estimated
% age of clinical transition, tau_hat_i (mean of 10 models). Every visit is
% placed on the disease-time axis
%       t*_ij = age at visit - tau_hat_i,
% so t* = 0 is the transition and t* < 0 the years before it.
cfg = load_cohort(COHORT, 'dataDir', C.dataDir, 'tauDir', C.tauDir);
D = cfg(TRANS).amyloid;                                % visits with baseline amyloid status
disp(D(1:5, {'SubjID', 'AGE_visit', 'Time_to_Event', MARKER, 'A_Status_Basal'}));

%% 2. Why a dedicated method? Three features of the data
% (a) Unequal number of visits per subject.
% (b) The subgroups occupy different regions of the axis: A+ subjects are
%     closer to the transition.
% (c) Visits of the same subject are not independent.
ok = ~isnan(D.(MARKER)) & ~isnan(D.Time_to_Event);
D  = D(ok, :);
g0 = D(~D.A_Status_Basal, :);   g1 = D(D.A_Status_Basal, :);
u0 = unique(g0.SubjID);  u1 = unique(g1.SubjID);
nv0 = arrayfun(@(s) sum(g0.SubjID == s), u0);
nv1 = arrayfun(@(s) sum(g1.SubjID == s), u1);
med0 = arrayfun(@(s) median(g0.Time_to_Event(g0.SubjID == s)), u0);
med1 = arrayfun(@(s) median(g1.Time_to_Event(g1.SubjID == s)), u1);
fprintf('A-: %d subjects, visits per subject %d-%d, median location %.1f years\n', ...
    numel(u0), min(nv0), max(nv0), median(med0));
fprintf('A+: %d subjects, visits per subject %d-%d, median location %.1f years\n', ...
    numel(u1), min(nv1), max(nv1), median(med1));
if DO_PLOTS
    figure('Color', 'w');
    histogram(med0, 'BinWidth', 1, 'FaceColor', [0 .6 0]); hold on;
    histogram(med1, 'BinWidth', 1, 'FaceColor', [.9 0 0]);
    xlabel('Median disease time of each subject (years)'); ylabel('Subjects');
    legend('A-', 'A+'); title('Location of each subgroup on the axis');
end

%% 3. Coverage range of each subgroup
% Each curve is drawn only where its subgroup has data: between the 10th
% and 90th percentiles of the subject medians.
tr0 = coverage_range(g0, RANGE_P, T_RANGE);
tr1 = coverage_range(g1, RANGE_P, T_RANGE);
fprintf('Coverage A-: [%.1f, %.1f] | A+: [%.1f, %.1f] years\n', tr0, tr1);

%% 4. Subject-weighted local-linear smoother (sw-LOESS)
% At each point t of the axis, every subject contributes one local weight,
% w_i(t) = max_j K((t*_ij - t)/bw), shared among his visits. A weighted
% line is fitted to the nearby visits and the curve is its value at t.
% No global shape is imposed, and a subject with many visits does not
% count more than one with few.
x0 = linspace(tr0(1), tr0(2), NQ);
x1 = linspace(tr1(1), tr1(2), NQ);
y0 = sw_loess_ll(g0.Time_to_Event, g0.(MARKER), g0.SubjID, u0, x0, BW);
y1 = sw_loess_ll(g1.Time_to_Event, g1.(MARKER), g1.SubjID, u1, x1, BW);
% For comparison, the usual global quadratic polynomial of each subgroup
p0 = polyfit(g0.Time_to_Event, g0.(MARKER), 2);
p1 = polyfit(g1.Time_to_Event, g1.(MARKER), 2);
if DO_PLOTS
    figure('Color', 'w'); hold on; box on; grid on;
    plot(g0.Time_to_Event, g0.(MARKER), '.', 'Color', [0 .6 0 .3]);
    plot(g1.Time_to_Event, g1.(MARKER), '.', 'Color', [.9 0 0 .3]);
    plot(x0, y0, '-', 'Color', [0 .5 0], 'LineWidth', 3);
    plot(x1, y1, '-', 'Color', [.8 0 0], 'LineWidth', 3);
    xx = linspace(T_RANGE(1), T_RANGE(2), NQ);
    plot(xx, polyval(p0, xx), '--', 'Color', [0 .5 0]);
    plot(xx, polyval(p1, xx), '--', 'Color', [.8 0 0]);
    xline(0, ':k');  xlim(T_RANGE);
    xlabel('Disease time t* (years)'); ylabel(MARKER);
    legend('A- visits', 'A+ visits', 'A- sw-LOESS', 'A+ sw-LOESS', ...
           'A- polynomial', 'A+ polynomial', 'Location', 'southwest');
    title('sw-LOESS on the coverage range vs. global polynomial');
end

%% 5. Test statistic on the common range
% T = mean |y1(t) - y0(t)| over the range covered by both subgroups.
t_lo = max(tr0(1), tr1(1));  t_hi = min(tr0(2), tr1(2));
xq = linspace(t_lo, t_hi, NQ);
uAll = [u0; u1];
Dall = [g0; g1];
curve = @(ids) sw_loess_ll(Dall.Time_to_Event, Dall.(MARKER), Dall.SubjID, ids, xq, BW);
T_obs = mean(abs(curve(u1) - curve(u0)), 'omitnan');
fprintf('Common range [%.1f, %.1f] years, T_obs = %.3f\n', t_lo, t_hi, T_obs);

%% 6. Subject-level permutation test
% The A+/A- labels are shuffled between SUBJECTS (all the visits of a
% subject keep the same label), and T is recomputed. The common range is
% kept fixed.   p = (1 + #{T_perm >= T_obs}) / (1 + B).
B = 199;  rng(C.SEED, 'twister');
T_perm = zeros(B, 1);
n1 = numel(u1);
for b = 1:B
    idx = randperm(numel(uAll));
    T_perm(b) = mean(abs(curve(uAll(idx(1:n1))) - curve(uAll(idx(n1+1:end)))), 'omitnan');
end
p_perm = (1 + sum(T_perm >= T_obs)) / (1 + B);
fprintf('Permutation p-value (%d permutations): %.3f\n', B, p_perm);

%% 7. Simulation-calibrated p-value
% The permutation test assumes that the labels are exchangeable. Here they
% are not: A+ subjects sit closer to the transition, and each tau_hat_i
% has an error that shifts all the visits of the subject. Together they
% can separate the two curves even when the true trajectory is common.
% The calibrated p-value simulates R null replicates with the same design
% (sizes, locations, visits, common curve, variances) and the
% transition-age error of the cohort (Table 1), and compares T_obs with
% the null distribution of T.
tmp = fullfile(tempdir, 'tutorial_framework');
bw_sensitivity(COHORT, 'estimator', 'linear', 'bwGrid', BW, 'analyses', {'A'}, ...
    'exportCalib', true, 'bCV', false, 'nPerm', B, 'outDir', tmp, ...
    'dataDir', C.dataDir, 'tauDir', C.tauDir);
K = readtable(fullfile(tmp, sprintf('calib_%s.csv', COHORT)), 'TextType', 'string');
TRS = {'CU2MCI', 'MCI2Dem'};
K = K(K.transition == TRS{TRANS} & K.var == MARKER, :);
writetable(K, fullfile(tmp, 'calib_one.csv'));
R = sim_framework('design', 'calibrated', 'calibFiles', {fullfile(tmp, 'calib_one.csv')}, ...
    'obsFiles', {fullfile(tmp, sprintf('bw_sensitivity_%s_ll.csv', COHORT))}, ...
    'tauTable', C.TAU_PAPER, 'nRepsCal', 100, 'nPerm', 99, 'bw', BW, ...
    'bPlots', false, 'outDir', fullfile(tmp, 'sim'));
r = R(strcmp({R.estimator}, sprintf('swll_bw%g', BW)));
fprintf(['Calibrated simulation: type I error of the permutation test = %.2f, ' ...
         'p_cal = %.3f\n'], r.reject, r.pCal);
% A difference is reported as robust when p < 0.05 and p_cal < 0.05.
