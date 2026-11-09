function R = bw_sensitivity(cohort, varargin)
%BW_SENSITIVITY  Sensitivity of the sw-LOESS kinematics to the bandwidth bw.
%
%   R = bw_sensitivity(cohort, Name, Value, ...)
%
%   Repeats, for every bw in a grid, the permutation tests of the
%   kinematics (Table 4 of the manuscript), for the three
%   stratifications:
%       'A'      A+ vs A-                  (row 2 of plot_tau_AV)
%       'V'      V+ vs V-                  (row 3 of plot_tau_AV)
%       'AxV_A+' V+ vs V- within A+        (row 1 of plot_tau_AV_interaction)
%       'AxV_A-' V+ vs V- within A-        (row 2 of plot_tau_AV_interaction)
%   Only sw-LOESS curves depend on bw. The marker that defines each
%   stratification (polynomial curve) is not included.
%
%   For each comparison, marker, transition and bw it stores:
%     - T_obs and the permutation p-value (Eqs. 7-8 of the manuscript);
%     - start/end values and slope of each subgroup curve (as in
%       Tables S6, S11 and S12);
%     - the minimum effective number of subjects along the common range,
%       n_eff(t) = (sum_i w_i(t))^2 / sum_i w_i(t)^2 (Kish);
%     - the leave-one-subject-out cross-validation error (CV-MSE).
%
%   Data are read with load_cohort (CSV files of data/). sw-LOESS and
%   the permutation test are those of plot_kinematics.m. The permutation test is vectorised: the subject-level terms
%   w_i(t) and ybar_i(t) do not depend on the labels, so they are
%   computed once and every permutation is a matrix product. The same
%   permutations (fixed seed) are reused for every bw, so differences
%   between bw are not Monte Carlo noise.
%
%   Name-value options
%     'bwGrid'     [2 3 4 6 8]   bandwidths (years)
%     'nPerm'      1000          permutations
%     'seed'       2026          RNG seed for the permutations
%     'bCV'        true          compute LOSO cross-validation
%     'invalidPerm' 'zero'       permutation with < 3 valid points:
%                                'zero' -> T_perm = 0 (as in the original
%                                scripts); 'drop' -> discarded
%     'outDir'     '.'
%     'bCheck'     true          check that the vectorised curves equal
%                                sw_loess_eval / sw_loess_ll_eval (max |diff| < 1e-9)
%     'estimator'  'const'       'const'  = sw-LOESS of the manuscript (local
%                                           constant, Eqs. 2-4);
%                                'linear' = subject-weighted local-linear
%                                           (sw_loess_ll_eval). Outputs get the
%                                           suffix _ll.
%     'exportCalib' false        also write calib_<cohort>.csv with the
%                                descriptors that sim_framework('design',
%                                'calibrated') needs (see export_calib).
%     'analyses'   {'A','V'}     comparisons: 'A', 'V', 'AxV_A+', 'AxV_A-'
%                                (paper A uses only 'A' and 'V')
%     'dataDir'    'data'        folder of the visit CSV files
%     'tauDir'     'data/tau'    folder of the transition-age CSV files
%     'wmhThreshold' 10          V+ threshold on WMH_bl (ml)
%
%   Outputs (in outDir)
%     bw_sensitivity_<cohort>.mat   table R with every result
%     bw_sensitivity_<cohort>.csv   same table
%     bw_pvalues_<cohort>.tex       p-values by bw (Supplementary table)
%     bw_cv_<cohort>.tex            CV-optimal bw and minimum n_eff
%     bw_summary_<cohort>.txt       stability of the four key results
%
%   The permutations use a fixed seed, so the p-values are reproduced
%   exactly with the same data, nPerm and seed.
%
%   Examples
%     R = bw_sensitivity('ADNI');
%     R = bw_sensitivity('NACC', 'bwGrid', [2 3 4 6 8], 'nPerm', 1000);
%
%   Authorship: Carlos Platero (UPM), with the collaboration of Claude
%   (Anthropic).

%% ---- options ------------------------------------------------------------
p = inputParser;
addRequired (p, 'cohort', @(x) any(strcmpi(x, {'ADNI', 'NACC'})));
addParameter(p, 'bwGrid', [2 3 4 6 8], @(x) isnumeric(x) && all(x > 0));
addParameter(p, 'nPerm', 1000, @(x) isscalar(x) && x >= 1);
addParameter(p, 'seed', 2026, @isscalar);
addParameter(p, 'bCV', true, @(x) isscalar(x));
addParameter(p, 'invalidPerm', 'zero', @(x) any(strcmpi(x, {'zero', 'drop'})));
addParameter(p, 'outDir', '.');
addParameter(p, 'bCheck', true, @(x) isscalar(x));
addParameter(p, 'estimator', 'const', @(x) any(strcmpi(x, {'const', 'linear'})));
addParameter(p, 'exportCalib', false, @(x) isscalar(x));
addParameter(p, 'analyses', {'A', 'V'});
addParameter(p, 'dataDir', 'data');
addParameter(p, 'tauDir', fullfile('data', 'tau'));
addParameter(p, 'wmhThreshold', 10);
parse(p, cohort, varargin{:});
o = p.Results;
o.cohort = upper(char(o.cohort));
o.invalidPerm = lower(char(o.invalidPerm));
o.bCV = logical(o.bCV);  o.bCheck = logical(o.bCheck);
o.estimator = lower(char(o.estimator));
o.exportCalib = logical(o.exportCalib);
sfx = ''; if strcmp(o.estimator, 'linear'), sfx = '_ll'; end
fprintf('[bw] Estimator: %s\n', o.estimator);
if ~isfolder(o.outDir), mkdir(o.outDir); end

%% ---- fixed parameters (identical to plot_tau_AV / _interaction) ---------
T_RANGE        = [-10, 5];
N_QUERY        = 150;
SUBGRP_RANGE_P = 10;

% sw-LOESS markers of each stratification (the defining marker is
% polynomial and does not depend on bw).
ANALYSES = struct( ...
    'name',  {'A',                'V',                'AxV_A+',            'AxV_A-'}, ...
    'label', {'A$+$ vs.\ A$-$',   'V$+$ vs.\ V$-$',   'V$+$ vs.\ V$-$, within A$+$', ...
              'V$+$ vs.\ V$-$, within A$-$'}, ...
    'vars',  {{'MEM','EF','WMH'}, {'MEM','EF','CENTILOIDS'}, ...
              {'MEM','EF','CENTILOIDS'}, {'MEM','EF','CENTILOIDS'}});

ANALYSES = ANALYSES(ismember({ANALYSES.name}, cellstr(o.analyses)));

%% ---- data (CSV files, same loader as plot_kinematics) -------------------
fprintf('\n[bw] Loading %s data...\n', o.cohort);
cfg = load_cohort(o.cohort, 'dataDir', o.dataDir, 'tauDir', o.tauDir, ...
                  'wmhThreshold', o.wmhThreshold);
trans_labels = {'CU2MCI', 'MCI2Dem'};

rows = {};  calib = {};
for g = 1:numel(cfg)
    joint = build_joint(cfg(g).select, o.wmhThreshold);   % PET+WMH subgroup (A x V)
    for a = 1:numel(ANALYSES)
        an = ANALYSES(a);
        switch an.name
            case 'A'
                d0 = cfg(g).amyloid(cfg(g).amyloid.A_Status_Basal == 0, :);
                d1 = cfg(g).amyloid(cfg(g).amyloid.A_Status_Basal == 1, :);
            case 'V'
                d0 = cfg(g).wmh(cfg(g).wmh.A_Status_Basal == 0, :);
                d1 = cfg(g).wmh(cfg(g).wmh.A_Status_Basal == 1, :);
            case 'AxV_A+'
                s  = joint(joint.A_Status == 1, :);
                d0 = s(s.V_Status == 0, :);  d1 = s(s.V_Status == 1, :);
            case 'AxV_A-'
                s  = joint(joint.A_Status == 0, :);
                d0 = s(s.V_Status == 0, :);  d1 = s(s.V_Status == 1, :);
        end
        fprintf('[bw] %s %s %-7s: %d vs %d subjects\n', o.cohort, trans_labels{g}, ...
            an.name, numel(unique(d0.SubjID)), numel(unique(d1.SubjID)));

        % Coverage ranges and common grid do not depend on bw.
        tr0 = subgrp_range(d0, SUBGRP_RANGE_P, T_RANGE);
        tr1 = subgrp_range(d1, SUBGRP_RANGE_P, T_RANGE);
        t_lo = max(tr0(1), tr1(1));  t_hi = min(tr0(2), tr1(2));

        % Same permutations for every marker and bw of this comparison.
        all_data = [d0; d1];
        u_all = unique(all_data.SubjID);
        n0 = numel(unique(d0.SubjID));  n1 = numel(unique(d1.SubjID));
        in0 = ismember(u_all, unique(d0.SubjID));
        rng(o.seed + 100*g + a, 'twister');
        P0 = false(numel(u_all), o.nPerm);  P1 = P0;
        for b = 1:o.nPerm
            idx = randperm(numel(u_all));
            P0(idx(1:n0), b) = true;
            P1(idx(n0+1:n0+n1), b) = true;
        end

        for v = 1:numel(an.vars)
            vn = an.vars{v};
            if o.exportCalib
                calib{end+1} = export_calib(d0, d1, vn, o.cohort, trans_labels{g}, an.name, T_RANGE); %#ok<AGROW>
            end
            for bw = o.bwGrid(:)'
                r = struct('cohort', string(o.cohort), 'transition', string(trans_labels{g}), ...
                    'analysis', string(an.name), 'var', string(vn), 'bw', bw, ...
                    'n0', n0, 'n1', n1, 'T_obs', NaN, 'p', NaN, ...
                    't_start0', NaN, 'y_start0', NaN, 't_end0', NaN, 'y_end0', NaN, 'slope0', NaN, ...
                    't_start1', NaN, 'y_start1', NaN, 't_end1', NaN, 'y_end1', NaN, 'slope1', NaN, ...
                    'neff_min0', NaN, 'neff_min1', NaN, 'cv_mse', NaN, 'cv_nsubj', NaN);

                % ---- permutation test (vectorised, Eqs. 7-8) ----
                if t_hi > t_lo
                    xq = linspace(t_lo, t_hi, N_QUERY);
                    Tm = make_terms(all_data, vn, u_all, xq, bw, o.estimator);
                    y0 = terms_curve(Tm, in0);
                    y1 = terms_curve(Tm, ~in0);
                    if o.bCheck && bw == o.bwGrid(1) && v == 1
                        check_against_original(d0, vn, xq, bw, y0, o.estimator);
                    end
                    ok = ~isnan(y0) & ~isnan(y1);
                    if sum(ok) >= 5
                        r.T_obs = mean(abs(y1(ok) - y0(ok)));
                        Y0p = terms_curve(Tm, P0);   % nPerm x nq
                        Y1p = terms_curve(Tm, P1);
                        okp = ~isnan(Y0p) & ~isnan(Y1p);
                        D = abs(Y1p - Y0p);  D(~okp) = 0;
                        nok = sum(okp, 2);
                        Tp = sum(D, 2) ./ max(nok, 1);
                        bad = nok < 3;
                        if strcmp(o.invalidPerm, 'zero')
                            Tp(bad) = 0;
                            r.p = (1 + sum(Tp >= r.T_obs)) / (1 + o.nPerm);
                        else
                            Tp = Tp(~bad);
                            r.p = (1 + sum(Tp >= r.T_obs)) / (1 + numel(Tp));
                        end
                    end
                    % Kish effective number of subjects on the common range
                    r.neff_min0 = min(kish(Tm.W(in0, :)));
                    r.neff_min1 = min(kish(Tm.W(~in0, :)));
                end

                % ---- plotted curves: start/end values and slopes ----
                [r.t_start0, r.y_start0, r.t_end0, r.y_end0, r.slope0] = ...
                    curve_metrics(d0, vn, bw, tr0, o.estimator, N_QUERY);
                [r.t_start1, r.y_start1, r.t_end1, r.y_end1, r.slope1] = ...
                    curve_metrics(d1, vn, bw, tr1, o.estimator, N_QUERY);

                % ---- leave-one-subject-out cross-validation ----
                if o.bCV
                    [e0, n0cv] = loso_cv(d0, vn, bw, T_RANGE, o.estimator);
                    [e1, n1cv] = loso_cv(d1, vn, bw, T_RANGE, o.estimator);
                    r.cv_nsubj = n0cv + n1cv;
                    if r.cv_nsubj > 0
                        r.cv_mse = (e0*n0cv + e1*n1cv) / r.cv_nsubj;
                    end
                end
                rows{end+1} = r; %#ok<AGROW>
            end
            fprintf('    %-10s p(bw) = %s\n', vn, ...
                strjoin(arrayfun(@(k) sprintf('%.3f', rows{end-numel(o.bwGrid)+k}.p), ...
                1:numel(o.bwGrid), 'UniformOutput', false), '  '));
        end
    end
end
R = struct2table([rows{:}]);

%% ---- outputs --------------------------------------------------------------
base = fullfile(o.outDir, ['bw_sensitivity_' o.cohort sfx]);
try
    save([base '.mat'], 'R', 'o');
catch err                                   % Octave cannot save table objects
    warning('bw_sensitivity:save', 'MAT file not written (%s); the CSV has the same table.', err.message);
end
writetable(R, [base '.csv']);
write_pvalue_table(R, o, ANALYSES, fullfile(o.outDir, ['bw_pvalues_' o.cohort sfx '.tex']));
if o.bCV
    write_cv_table(R, o, ANALYSES, fullfile(o.outDir, ['bw_cv_' o.cohort sfx '.tex']));
end
write_summary(R, o, ANALYSES, fullfile(o.outDir, ['bw_summary_' o.cohort sfx '.txt']));
if o.exportCalib
    C = struct2table([calib{:}]);
    writetable(C, fullfile(o.outDir, ['calib_' o.cohort '.csv']));
    fprintf('[bw] Calibration descriptors saved to calib_%s.csv\n', o.cohort);
end
fprintf('\n[bw] Results saved with prefix %s\n', base);
end


%% =========================================================================
%  PET+WMH subgroup (same rule as plot_tau_AV_interaction)
%% =========================================================================
function joint = build_joint(sel, wmhThreshold)
% A_bl already holds the cohort rule (ADNI: PET status; NACC: CL > 20).
joint = sel(~isnan(sel.A_bl) & ~isnan(sel.WMH_bl), :);
joint.A_Status = joint.A_bl > 0;
joint.V_Status = joint.WMH_bl > wmhThreshold;
end


%% =========================================================================
%  Vectorised sw-LOESS (Eqs. 2-4)
%% =========================================================================
function n = kish(Wsub)
s1 = sum(Wsub, 1);  s2 = sum(Wsub.^2, 1);
n = (s1.^2) ./ s2;
n(s2 == 0) = 0;
end


function Tm = make_terms(data, vn, u_all, xq, bw, est)
% Per-subject terms on the query grid. A subset curve is a function of
% sums of these terms over the subjects of the subset (terms_curve), so
% every permutation is a matrix product.
%   'const'  (sw-LOESS):          W = w_i(t), WY = w_i(t) * ybar_i(t)
%   'linear' (sw local-linear):   weighted moments of the local line with
%            visit weights a_ij = w_i(t) k_ij(t) / sum_j k_ij(t)
n_s = numel(u_all);  n_q = numel(xq);
Tm.est = est;
Tm.W = zeros(n_s, n_q);  Tm.hasD = false(n_s, 1);
if strcmp(est, 'linear')
    Tm.S1 = Tm.W; Tm.S2 = Tm.W; Tm.T0 = Tm.W; Tm.T1 = Tm.W;
else
    Tm.WY = Tm.W;
end
valid = ~isnan(data.(vn)) & ~isnan(data.Time_to_Event);
t_all = data.Time_to_Event(valid);  y_all = data.(vn)(valid);  id_all = data.SubjID(valid);
[tf, loc] = ismember(id_all, u_all);
t_all = t_all(tf);  y_all = y_all(tf);  loc = loc(tf);
for si = unique(loc)'
    m = loc == si;
    dt = t_all(m) - xq;                               % n_vis x n_q
    K = exp(-0.5 * (dt / bw).^2);
    w = max(K, [], 1);
    w(w < 1e-12) = 0;
    Tm.hasD(si) = true;
    if strcmp(est, 'linear')
        a = K .* (w ./ max(sum(K, 1), realmin));
        ys = y_all(m);
        Tm.W(si, :)  = sum(a, 1);           Tm.S1(si, :) = sum(a .* dt, 1);
        Tm.S2(si, :) = sum(a .* dt.^2, 1);  Tm.T0(si, :) = sum(a .* ys, 1);
        Tm.T1(si, :) = sum(a .* dt .* ys, 1);
    else
        ybar = sum(K .* y_all(m), 1) ./ sum(K, 1);
        Tm.W(si, :)  = w;
        wy = w .* ybar;  wy(w == 0) = 0;
        Tm.WY(si, :) = wy;
    end
end
end


function Y = terms_curve(Tm, M)
% M: n_s x k (logical or counts), one column per subset. Y: k x n_q.
% A subset with < 3 subjects with data gives NaN (as eval_curve).
Mf = double(M)';
S0 = Mf * Tm.W;
if strcmp(Tm.est, 'linear')
    S1 = Mf * Tm.S1;  S2 = Mf * Tm.S2;  T0 = Mf * Tm.T0;  T1 = Mf * Tm.T1;
    den = S0 .* S2 - S1.^2;
    Y = (S2 .* T0 - S1 .* T1) ./ den;
    Y(den <= 1e-10 * max(S0.^2, realmin)) = NaN;
else
    Y = (Mf * Tm.WY) ./ S0;
end
Y(S0 <= 1e-12) = NaN;
few = sum((M ~= 0) & Tm.hasD, 1)' < 3;
Y(few, :) = NaN;
end


function c = export_calib(d0, d1, vn, cohort, trans, analysis, T_RANGE)
% Descriptors of one comparison for sim_framework('design','calibrated'):
% sizes, location of each subgroup on the disease-time axis (subject
% medians), visit pattern, and the global quadratic curve of the marker
% with its between- and within-subject SDs (residuals of that curve).
c = struct('cohort', string(cohort), 'transition', string(trans), ...
    'analysis', string(analysis), 'var', string(vn));
D = {d0, d1};
for g = 0:1
    d = D{g+1};
    ok = ~isnan(d.(vn)) & ~isnan(d.Time_to_Event);
    ids = d.SubjID(ok);  t = d.Time_to_Event(ok);
    u = unique(ids);
    med = zeros(numel(u), 1);  nv = med;
    for k = 1:numel(u)
        m = ids == u(k);  med(k) = median(t(m));  nv(k) = sum(m);
    end
    c.(sprintf('n%d', g))       = numel(u);
    c.(sprintf('mu%d', g))      = mean(med);
    c.(sprintf('sd%d', g))      = std(med);
    c.(sprintf('pSingle%d', g)) = mean(nv == 1);
    c.(sprintf('nvMulti%d', g)) = mean(nv(nv > 1));
    c.(sprintf('gap%d', g))     = NaN;
    gaps = [];
    for k = 1:numel(u)
        tk = sort(t(ids == u(k)));
        if numel(tk) > 1, gaps = [gaps; diff(tk)]; end %#ok<AGROW>
    end
    if ~isempty(gaps), c.(sprintf('gap%d', g)) = median(gaps); end
end
% Global quadratic curve and variance components (both subgroups pooled).
d = [d0; d1];
ok = ~isnan(d.(vn)) & ~isnan(d.Time_to_Event) & ...
     d.Time_to_Event >= T_RANGE(1) & d.Time_to_Event <= T_RANGE(2);
t = d.Time_to_Event(ok);  y = d.(vn)(ok);  ids = d.SubjID(ok);
pp = polyfit(t, y, 2);
c.c2 = pp(1);  c.c1 = pp(2);  c.c0 = pp(3);
res = y - polyval(pp, t);
u = unique(ids);
mres = zeros(numel(u), 1);  nv = mres;  ssw = 0;  dfw = 0;
for k = 1:numel(u)
    m = ids == u(k);
    mres(k) = mean(res(m));  nv(k) = sum(m);
    if nv(k) > 1
        ssw = ssw + sum((res(m) - mres(k)).^2);  dfw = dfw + nv(k) - 1;
    end
end
sigE = sqrt(ssw / max(dfw, 1));
sigB = sqrt(max(var(mres) - mean(sigE^2 ./ nv), 0));
c.sigB = sigB;  c.sigE = sigE;
end


function y_hat = sw_loess_ll_eval(t, y, ids, u_ids, x_plot, bw)
% Subject-weighted local-linear smoother (reference implementation; same
% as the stand-alone sw_loess_ll_eval.m).
n_q = numel(x_plot);
y_hat = nan(1, n_q);
S0 = zeros(1, n_q); S1 = S0; S2 = S0; T0 = S0; T1 = S0;
x_plot = x_plot(:)';
for si = 1:numel(u_ids)
    m = ids == u_ids(si);
    if ~any(m), continue; end
    dt = t(m) - x_plot;
    K  = exp(-0.5 * (dt / bw).^2);
    w  = max(K, [], 1);  w(w < 1e-12) = 0;
    a  = K .* (w ./ max(sum(K, 1), realmin));
    ys = y(m);
    S0 = S0 + sum(a, 1);           S1 = S1 + sum(a .* dt, 1);
    S2 = S2 + sum(a .* dt.^2, 1);  T0 = T0 + sum(a .* ys, 1);
    T1 = T1 + sum(a .* dt .* ys, 1);
end
den = S0 .* S2 - S1.^2;
ok = den > 1e-10 * max(S0.^2, realmin) & S0 > 1e-12;
y_hat(ok) = (S2(ok) .* T0(ok) - S1(ok) .* T1(ok)) ./ den(ok);
end


function check_against_original(d, vn, xq, bw, y_vec, est)
valid = ~isnan(d.(vn)) & ~isnan(d.Time_to_Event);
u = unique(d.SubjID(valid));
if numel(u) < 3, return; end
y_ref = eval_curve_est(d.Time_to_Event(valid), d.(vn)(valid), d.SubjID(valid), u, xq, bw, est);
dd = abs(y_ref - y_vec);  dif = max([0, dd(~isnan(dd))]);
if dif > 1e-9
    error('bw_sensitivity:check', 'Vectorised curve differs from the reference smoother (%.2e).', dif);
end
fprintf('    [check] vectorised %s smoother = reference (max diff %.1e)\n', est, dif);
end


function y = eval_curve_est(t, yv, ids, u, x, bw, est)
if strcmp(est, 'linear')
    y = sw_loess_ll_eval(t, yv, ids, u, x, bw);
else
    y = sw_loess_eval(t, yv, ids, u, x, bw);
end
end


%% =========================================================================
%  Plotted-curve metrics (as plot_trend_loess_sw_out + make_metric)
%% =========================================================================
function [ts, ys, te, ye, sl] = curve_metrics(data, vn, bw, tr, est, n_q)
ts = NaN; ys = NaN; te = NaN; ye = NaN; sl = NaN;
valid = ~isnan(data.(vn)) & ~isnan(data.Time_to_Event);
if sum(valid) < 10, return; end
t_all = data.Time_to_Event(valid);  y_all = data.(vn)(valid);  id_all = data.SubjID(valid);
u = unique(id_all);
if numel(u) < 5, return; end
t_lo = max(tr(1), prctile(t_all, 5));
t_hi = min(tr(2), prctile(t_all, 95));
if t_hi <= t_lo, return; end
x = linspace(t_lo, t_hi, n_q);
y = eval_curve_est(t_all, y_all, id_all, u, x, bw, est);
if all(isnan(y)), return; end
ts = x(1); te = x(end); ys = y(1); ye = y(end);
sl = (ye - ys) / (te - ts);
end


%% =========================================================================
%  Leave-one-subject-out cross-validation
%% =========================================================================
function [cv, n_used] = loso_cv(data, vn, bw, T_RANGE, est)
% For each subject i, predict its visits inside T_RANGE with the curve
% of the other subjects of the same subgroup. Error of subject i: mean
% squared error of its visits. CV: mean over subjects (one weight per
% subject, consistent with sw-LOESS). The leave-one-out curves of all
% subjects are obtained at once with the subset matrix ~eye(n).
cv = NaN; n_used = 0;
valid = ~isnan(data.(vn)) & ~isnan(data.Time_to_Event) & ...
        data.Time_to_Event >= T_RANGE(1) & data.Time_to_Event <= T_RANGE(2);
if sum(valid) < 10, return; end
t_all = data.Time_to_Event(valid);  y_all = data.(vn)(valid);  id_all = data.SubjID(valid);
u = unique(id_all);
if numel(u) < 5, return; end
xg = linspace(T_RANGE(1), T_RANGE(2), 301);
Tm = make_terms(table(id_all, t_all, y_all, ...
          'VariableNames', {'SubjID', 'Time_to_Event', vn}), vn, u, xg, bw, est);
Yloo = terms_curve(Tm, ~eye(numel(u)));       % n x numel(xg)
err = nan(numel(u), 1);
[~, loc] = ismember(id_all, u);
for si = 1:numel(u)
    m = loc == si;
    yp = interp1(xg, Yloo(si, :), t_all(m), 'linear');
    e = (y_all(m) - yp).^2;
    if any(~isnan(e)), err(si) = mean(e, 'omitnan'); end
end
n_used = sum(~isnan(err));
if n_used > 0, cv = mean(err, 'omitnan'); end
end


%% =========================================================================
%  LaTeX tables (Supplementary Material style)
%% =========================================================================
function write_pvalue_table(R, o, AN, file)
fid = fopen(file, 'w', 'n', 'UTF-8');  c = onCleanup(@() fclose(fid));
wl = @(s) fprintf(fid, '%s\n', s);
bws = o.bwGrid(:)';  nb = numel(bws);
varlab = containers.Map({'MEM','EF','CENTILOIDS','WMH'}, {'MEM','EF','Centiloids','WMH'});
wl(['% Auto-generated by bw_sensitivity.m (' o.cohort ', nPerm = ' num2str(o.nPerm) ...
    ', seed = ' num2str(o.seed) ', invalidPerm = ' o.invalidPerm ')']);
wl('\begin{table}[htbp]');
wl('\centering');
wl(['\caption{Sensitivity of the permutation tests to the bandwidth $bw$ of the \textit{sw}-LOESS smoother, ' o.cohort ...
    '. Each cell is the $p$-value of the comparison between the two subgroups. The permutations are the same for every $bw$. ' ...
    'The marker that defines each stratification is not included. ' ...
    'Bold: $p<0.05$. n.e.: not estimable.}']);
wl(['\label{tab:S_bw_' lower(o.cohort) '}']);
wl('\scriptsize');
wl('\setlength{\tabcolsep}{3pt}');
wl(['\begin{tabular}{@{}ll' repmat('c', 1, 2*nb) '@{}}']);
wl('\toprule');
wl(['& & \multicolumn{' num2str(nb) '}{c}{CU$\to$MCI} & \multicolumn{' num2str(nb) '}{c}{MCI$\to$Dem.} \\']);
wl(sprintf('\\cmidrule(lr){3-%d}\\cmidrule(lr){%d-%d}', 2+nb, 3+nb, 2+2*nb));
hdr = strjoin(arrayfun(@(b) sprintf('$bw{=}%g$', b), bws, 'UniformOutput', false), ' & ');
wl(['\textbf{Comparison} & \textbf{Marker} & ' hdr ' & ' hdr ' \\']);
wl('\midrule');
for a = 1:numel(AN)
    for v = 1:numel(AN(a).vars)
        vn = AN(a).vars{v};
        if v == 1, lead = AN(a).label; else, lead = ''; end
        cells = {};
        for trc = {'CU2MCI', 'MCI2Dem'}, tr = trc{1};
            for b = bws
                k = R.analysis == AN(a).name & R.var == vn & R.transition == tr & R.bw == b;
                cells{end+1} = fmt_p(R.p(k)); %#ok<AGROW>
            end
        end
        wl([lead ' & ' varlab(vn) ' & ' strjoin(cells, ' & ') ' \\']);
    end
    if a < numel(AN), wl('\midrule'); end
end
wl('\bottomrule');
wl('\end{tabular}');
wl('\end{table}');
end


function write_cv_table(R, o, AN, file)
fid = fopen(file, 'w', 'n', 'UTF-8');  c = onCleanup(@() fclose(fid));
wl = @(s) fprintf(fid, '%s\n', s);
varlab = containers.Map({'MEM','EF','CENTILOIDS','WMH'}, {'MEM','EF','Centiloids','WMH'});
wl(['% Auto-generated by bw_sensitivity.m (' o.cohort ')']);
wl('\begin{table}[htbp]');
wl('\centering');
wl(['\caption{Bandwidth selected by leave-one-subject-out cross-validation, ' o.cohort ...
    '. $bw_{\rm CV}$: grid value with the smallest mean squared prediction error. ' ...
    '$n_{\rm eff}$: minimum effective number of subjects along the common range, in the smaller subgroup, ' ...
    'with $bw=6$ and with $bw_{\rm CV}$. Grid: ' strjoin(arrayfun(@(b) sprintf('%g', b), o.bwGrid, 'UniformOutput', false), ', ') ' years.}']);
wl(['\label{tab:S_bwcv_' lower(o.cohort) '}']);
wl('\scriptsize');
wl('\begin{tabular}{@{}llcccccc@{}}');
wl('\toprule');
wl('& & \multicolumn{3}{c}{CU$\to$MCI} & \multicolumn{3}{c}{MCI$\to$Dem.} \\');
wl('\cmidrule(lr){3-5}\cmidrule(lr){6-8}');
wl(['\textbf{Comparison} & \textbf{Marker} & $bw_{\rm CV}$ & $n_{\rm eff}(6)$ & $n_{\rm eff}(bw_{\rm CV})$ ' ...
    '& $bw_{\rm CV}$ & $n_{\rm eff}(6)$ & $n_{\rm eff}(bw_{\rm CV})$ \\']);
wl('\midrule');
for a = 1:numel(AN)
    for v = 1:numel(AN(a).vars)
        vn = AN(a).vars{v};
        if v == 1, lead = AN(a).label; else, lead = ''; end
        cells = {};
        for trc = {'CU2MCI', 'MCI2Dem'}, tr = trc{1};
            k = R.analysis == AN(a).name & R.var == vn & R.transition == tr;
            Rk = R(k, :);
            [~, i] = min(Rk.cv_mse);
            if isempty(i) || all(isnan(Rk.cv_mse))
                cells = [cells, {'n.e.', 'n.e.', 'n.e.'}]; %#ok<AGROW>
                continue;
            end
            neff = min(Rk.neff_min0, Rk.neff_min1);
            n6 = neff(Rk.bw == 6);
            if isempty(n6), s6 = '--'; else, s6 = sprintf('%.0f', n6); end
            cells = [cells, {sprintf('%g', Rk.bw(i)), s6, sprintf('%.0f', neff(i))}]; %#ok<AGROW>
        end
        wl([lead ' & ' varlab(vn) ' & ' strjoin(cells, ' & ') ' \\']);
    end
    if a < numel(AN), wl('\midrule'); end
end
wl('\bottomrule');
wl('\end{tabular}');
wl('\end{table}');
end


function s = fmt_p(p)
if isempty(p) || isnan(p), s = 'n.e.'; return; end
if p < 0.001, s = '$<$0.001'; else, s = sprintf('%.3f', p); end
if p < 0.05, s = ['\textbf{' s '}']; end
end


%% =========================================================================
%  Stability of the four key results
%% =========================================================================
function write_summary(R, o, AN, file)
key = { 'A',      'MEM', 'significant', 'MEM, A+ vs A-';
        'V',      'EF',  'significant', 'EF, V+ vs V-';
        'AxV_A+', 'EF',  'significant', 'EF, V+ vs V- within A+';
        'AxV_A-', 'EF',  'null',        'EF, V+ vs V- within A- (no effect)' };
key = key(ismember(key(:, 1), {AN.name}), :);
fid = fopen(file, 'w');  c = onCleanup(@() fclose(fid));
out = @(varargin) emit(fid, varargin{:});
out('Bandwidth sensitivity, %s (nPerm = %d, seed = %d)\n', o.cohort, o.nPerm, o.seed);
out('Grid: %s\n\n', mat2str(o.bwGrid));
for trc = {'CU2MCI', 'MCI2Dem'}, tr = trc{1};
    out('%s\n', tr);
    for k = 1:size(key, 1)
        m = R.analysis == key{k,1} & R.var == key{k,2} & R.transition == tr;
        pv = R.p(m);  bw = R.bw(m);
        sig = pv < 0.05;
        if strcmp(key{k,3}, 'significant'), expected = sig; else, expected = ~sig; end
        if all(expected),      status = 'STABLE';
        elseif ~any(expected), status = 'NOT AS EXPECTED AT ANY bw';
        else,                  status = 'DEPENDS ON bw';
        end
        out('  %-38s %-26s p = %s (bw = %s)\n', key{k,4}, status, ...
            mat2str(round(1000 * pv') / 1000), mat2str(bw'));
    end
    out('\n');
end
end


function emit(fid, fmt, varargin)
% Write the same line to the console and to the summary file.
fprintf(1, fmt, varargin{:});
fprintf(fid, fmt, varargin{:});
end


%% =========================================================================
%  COPIED FROM plot_kinematics.m WITHOUT CHANGES IN LOGIC
%  (subgrp_range, sw_loess_eval)
%% =========================================================================
function tr = subgrp_range(data, range_p, t_range)
    valid = ~isnan(data.Time_to_Event);
    if sum(valid) < 5, tr = t_range; return; end
    u_ids  = unique(data.SubjID(valid));
    t_subj = arrayfun(@(id) median(data.Time_to_Event( ...
                 data.SubjID==id & valid)), u_ids);
    t_subj = t_subj(~isnan(t_subj));
    if numel(t_subj) < 3, tr = t_range; return; end
    t_lo = max(prctile(t_subj, range_p),       t_range(1));
    t_hi = min(prctile(t_subj, 100-range_p),   t_range(2));
    if t_hi <= t_lo, tr = t_range; return; end
    if (t_hi-t_lo) >= 0.90*(t_range(2)-t_range(1))
        tr = t_range;
    else
        tr = [t_lo, t_hi];
    end
end


function y_hat = sw_loess_eval(t, y, ids, u_ids, x_plot, bw)
    n_q = numel(x_plot);
    y_hat = nan(1, n_q);
    n_s = numel(u_ids);
    t_cell = cell(n_s,1); y_cell = cell(n_s,1);
    for si = 1:n_s
        m = ids == u_ids(si);
        t_cell{si} = t(m); y_cell{si} = y(m);
    end
    for qi = 1:n_q
        t0 = x_plot(qi); num = 0; denom = 0;
        for si = 1:n_s
            k_s = exp(-0.5*((t_cell{si}-t0)/bw).^2);
            if isempty(k_s), continue; end
            w_s = max(k_s);
            if w_s < 1e-12, continue; end
            ybar_s = sum(k_s.*y_cell{si}) / sum(k_s);
            num   = num   + w_s*ybar_s;
            denom = denom + w_s;
        end
        if denom > 1e-12, y_hat(qi) = num/denom; end
    end
end
