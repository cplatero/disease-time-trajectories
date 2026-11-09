function R = sim_framework(varargin)
% SIM_FRAMEWORK  Simulation study of the subject-weighted trajectory
% comparison framework (sw-LOESS + subject-level permutation test).
%
% Synthetic subjects are generated with a KNOWN population trajectory per
% subgroup. The framework, and two alternatives, are applied to each
% replicate, and their behaviour is measured:
%   * type I error of the permutation test (scenario 'null');
%   * power (scenarios 'level' and 'diverge');
%   * bias / RMSE of the estimated curves and of their difference;
%   * pointwise coverage of the subject-bootstrap 95% bands (sw-LOESS).
%
% The data reproduce the features of ADNI/NACC that motivate the framework:
%   1. unequal number of visits per subject, optionally informative (subjects
%      with better scores stay longer in the cohort and have more visits);
%   2. subgroups located in different regions of the disease-time axis
%      (reference subgroup far before the event, positive subgroup near it);
%   3. a fraction of subjects with a single measurement (NACC imaging);
%   4. error in the estimated transition age (shift of the whole subject).
%
% Estimators compared (all on the same data and the same permutations):
%   sw_bwX : sw-LOESS (Eqs. 2-4 of the manuscript) with bandwidth X years;
%   nw     : classical visit-weighted Nadaraya-Watson, bandwidth = 'bw';
%   poly   : global quadratic polynomial per subgroup (all visits in the
%            window), as in the usual natural-history plots.
% The test statistic is the same for all of them (Eq. 7): mean absolute
% difference between the two curves on the common coverage range, which
% is fixed across permutations. p = (1 + #{T_perm >= T_obs}) / (1 + nPerm).
%
% USAGE
%   R = sim_framework();                         % full study (long)
%   R = sim_framework('design', 'quick');        % smoke test, < 1 min
%   R = sim_framework('design', 'typeI', 'nRepsNull', 2000);
%   R = sim_framework('design', 'power', 'sigTau', 1.5, 'pSingle', 0.85);
%   R = sim_framework('design', 'calibrated', ...
%         'calibFiles', {'calib_ADNI.csv', 'calib_NACC.csv'}, ...
%         'obsFiles',   {'bw_sensitivity_ADNI.csv', 'bw_sensitivity_ADNI_ll.csv', ...
%                        'bw_sensitivity_NACC.csv', 'bw_sensitivity_NACC_ll.csv'});
%
% CALIBRATED DESIGN
%   One null scenario per real comparison and marker. The calib_*.csv files
%   (bw_sensitivity(..., 'exportCalib', true)) give, per subgroup, the
%   number of subjects, the location of the subject medians on the axis,
%   the visit pattern, and the global quadratic curve of the marker with
%   its between/within-subject SDs. Both subgroups share that curve (null
%   hypothesis). The transition-age error is a normal mixture per cohort
%   and transition ('tauTable'), fitted to the robust and classical SDs of
%   tau_hat - clinical onset in converters. For each comparison the run
%   gives the type I error and the null distribution of T. If obsFiles
%   (bw_sensitivity outputs, '_ll' = local-linear) are given, it also
%   reports the observed T at bw = 'bw' and a calibrated p-value,
%   pCal = (1 + #{T_null >= T_obs}) / (1 + nRepsCal), which accounts for
%   the design and the transition-age error.
%
% NAME-VALUE OPTIONS (defaults in brackets)
%   'design'     'all' | 'typeI' | 'power' | 'quick'           ['all']
%   'nGrid'      subjects per subgroup                         [25 75 300]
%   'effGrid'    effect sizes, in SD of the marker            [0.1 0.2 0.3 0.5]
%   'sigTauGrid' SD (years) of the transition-age error,
%                type I design only                           [0 1.5]
%   'pSingleGrid' fraction of single-visit subjects,
%                type I design only                           [0.15 0.85]
%   'sigTau'     transition-age error SD, power design        [1]
%   'pSingle'    single-visit fraction, power design          [0.15]
%   'bw'         main bandwidth (years), also used by 'nw'    [3]
%   'bwGrid'     sw-LOESS bandwidths evaluated                [2 3 4 6]
%   'nRepsNull'  replicates per type I cell                   [1000]
%   'nRepsAlt'   replicates per power cell                    [500]
%   'nPerm'      permutations per replicate                   [500]
%   'nBoot'      bootstrap resamples for the bands (0 = off)  [200]
%   'seed'       random seed                                  [2026]
%   'outDir'     output folder                                ['sim_out']
%   'bPlots'     draw figures                                 [true]
%   'gen'        struct overriding generative parameters (see default_gen)
%   'tauModel'   'normal' | 'mix' (typeI/power designs)            ['normal']
%   'tauCore', 'tauTail', 'pTail'  mixture for tauModel = 'mix'     [1 6 0.15]
%   'calibFiles' / 'obsFiles' / 'nRepsCal' / 'tauTable'   calibrated design
%   'calRobust'  true: centres not clipped and null curve flat 2 y beyond
%                the window (sensitivity check)                     [false]
%   'calNoTau'   true: calibrated replicates without transition-age
%                error (delta_i = 0); the true spread of the locations
%                is the same as in the main run (sensitivity check) [false]
%
% OUTPUTS (in outDir)
%   sim_results.mat / sim_results.csv : one row per cell and estimator
%   sim_typeI.tex, sim_power.tex, sim_bias.tex, sim_coverage.tex
%   sim_power.png, sim_typeI.png (if bPlots)
%
% Calibration. The defaults approximate a standardised cognitive composite
% (ADNI_MEM-like scale). sigTau should be set to the typical SD of the
% transition age across the 10 Leaspy models; it can be computed as
% mean(std(age_onset(:, 2:end), 0, 2, 'omitnan'), 'omitnan').
%
% Authorship: Carlos Platero (UPM), with the collaboration of Claude
% (Anthropic). Runs in MATLAB R2019b+ and GNU Octave 8.

%% ------------------------------------------------------------------ options
p = inputParser;
addParameter(p, 'design', 'all');
addParameter(p, 'nGrid', [25 75 300]);
addParameter(p, 'effGrid', [0.1 0.2 0.3 0.5]);
addParameter(p, 'sigTauGrid', [0 1.5]);
addParameter(p, 'pSingleGrid', [0.15 0.85]);
addParameter(p, 'sigTau', 1);
addParameter(p, 'pSingle', 0.15);
addParameter(p, 'bw', 3);
addParameter(p, 'bwGrid', [2 3 4 6]);
addParameter(p, 'nRepsNull', 1000);
addParameter(p, 'nRepsAlt', 500);
addParameter(p, 'nPerm', 500);
addParameter(p, 'nBoot', 200);
addParameter(p, 'seed', 2026);
addParameter(p, 'outDir', 'sim_out');
addParameter(p, 'bPlots', true);
addParameter(p, 'gen', struct());
addParameter(p, 'nStrata', 5);
addParameter(p, 'tauModel', 'normal');
addParameter(p, 'tauCore', 1);
addParameter(p, 'tauTail', 6);
addParameter(p, 'pTail', 0.15);
addParameter(p, 'calibFiles', {});
addParameter(p, 'obsFiles', {});
addParameter(p, 'nRepsCal', 500);
addParameter(p, 'tauTable', default_tau_table());
addParameter(p, 'calRobust', false);   % calibrated: no clipping, flat null beyond window
addParameter(p, 'calNoTau', false);    % calibrated: delta_i = 0 (design bias only)
parse(p, varargin{:});
o = p.Results;
o.design = lower(char(o.design));
o.bwGrid = unique([o.bwGrid(:)' o.bw]);

if strcmp(o.design, 'calibrated'), o.bwGrid = o.bw; end
if strcmp(o.design, 'quick')
    o.nGrid = [40 150]; o.effGrid = [0.3]; o.sigTauGrid = 1; o.pSingleGrid = 0.15;
    o.nRepsNull = 40; o.nRepsAlt = 20; o.nPerm = 199; o.nBoot = 50;
end

G = default_gen();
fn = fieldnames(o.gen);
for k = 1:numel(fn), G.(fn{k}) = o.gen.(fn{k}); end

if ~exist(o.outDir, 'dir'), mkdir(o.outDir); end
rng(o.seed, 'twister');

est = estimator_names(o);
fprintf('[sim] Estimators: %s\n', strjoin(est, ', '));

%% ------------------------------------------------------------------ cells
cells = struct('design', {}, 'scenario', {}, 'n', {}, 'eff', {}, ...
               'sigTau', {}, 'pSingle', {}, 'nReps', {}, 'cal', {}, 'tau', {});
tauDef = struct('model', o.tauModel, 'core', o.tauCore, 'tail', o.tauTail, 'pTail', o.pTail);
if any(strcmp(o.design, {'all', 'typei', 'quick'}))
    for n = o.nGrid
        for st = o.sigTauGrid
            for ps = o.pSingleGrid
                tz = tauDef; if strcmp(tz.model, 'normal'), tz.core = st; end
                cells(end+1) = struct('design', 'typeI', 'scenario', 'null', ...
                    'n', n, 'eff', 0, 'sigTau', st, 'pSingle', ps, ...
                    'nReps', o.nRepsNull, 'cal', [], 'tau', tz); %#ok<AGROW>
            end
        end
    end
end
if any(strcmp(o.design, {'all', 'power', 'quick'}))
    scen = {'level', 'diverge'};
    for s = 1:numel(scen)
        for n = o.nGrid
            for e = o.effGrid
                tz = tauDef; if strcmp(tz.model, 'normal'), tz.core = o.sigTau; end
                cells(end+1) = struct('design', 'power', 'scenario', scen{s}, ...
                    'n', n, 'eff', e, 'sigTau', o.sigTau, 'pSingle', o.pSingle, ...
                    'nReps', o.nRepsAlt, 'cal', [], 'tau', tz); %#ok<AGROW>
            end
        end
    end
end
if strcmp(o.design, 'calibrated')
    [cells, OBS] = calibrated_cells(o);
else
    OBS = [];
end
fprintf('[sim] %d cells, %d replicates in total\n', numel(cells), sum([cells.nReps]));

%% ------------------------------------------------------------------ run
rows = {};
t0 = tic;
for c = 1:numel(cells)
    C = cells(c);
    nr = C.nReps; ne = numel(est);
    pv   = nan(nr, ne);  ps   = nan(nr, ne);  bias = nan(nr, ne);  rmse = nan(nr, ne);
    dbia = nan(nr, ne);  drms = nan(nr, ne);  cov  = nan(nr, ne);  Tn = nan(nr, ne);
    doBoot = o.nBoot > 0 && strcmp(C.design, 'power');
    for r = 1:nr
        D = generate_data(C, G);
        M = analyse(D, G, o, est, doBoot);
        pv(r, :) = M.p(1, :);  ps(r, :) = M.p(2, :);
        bias(r, :) = M.bias;  rmse(r, :) = M.rmse;
        dbia(r, :) = M.dbias; drms(r, :) = M.drmse; cov(r, :) = M.cover; Tn(r, :) = M.T;
    end
    for k = 1:ne
        nv  = sum(~isnan(pv(:, k)));
        rej = sum(pv(:, k) < 0.05) / max(nv, 1);
        rjs = sum(ps(:, k) < 0.05) / max(nv, 1);
        tk = Tn(~isnan(Tn(:, k)), k);
        Tobs = lookup_obs(OBS, C, est{k}, o.bw);
        pCal = NaN;
        if ~isnan(Tobs) && ~isempty(tk), pCal = (1 + sum(tk >= Tobs)) / (1 + numel(tk)); end
        t95 = NaN; if ~isempty(tk), t95 = prctile(tk, 95); end
        rows(end+1, :) = {C.design, C.scenario, C.n, C.eff, C.sigTau, C.pSingle, ...
            est{k}, rej, sqrt(rej*(1-rej)/max(nv,1)), rjs, nv, ...
            mean(bias(:,k), 'omitnan'), mean(rmse(:,k), 'omitnan'), ...
            mean(dbia(:,k), 'omitnan'), mean(drms(:,k), 'omitnan'), ...
            mean(cov(:,k), 'omitnan'), mean(tk), t95, Tobs, pCal}; %#ok<AGROW>
    end
    rj = sum(pv < 0.05, 1) ./ max(sum(~isnan(pv), 1), 1);
    rs = sum(ps < 0.05, 1) ./ max(sum(~isnan(ps), 1), 1);
    txt = '';
    for k = 1:ne, txt = [txt sprintf('%s=%.2f/%.2f ', est{k}, rj(k), rs(k))]; end %#ok<AGROW>
    fprintf('[sim] cell %3d/%d  %-6s %-8s n=%4d eff=%.2f sigTau=%.1f pSingle=%.2f  | %s (%.0f s)\n', ...
        c, numel(cells), C.design, C.scenario, C.n, C.eff, C.sigTau, C.pSingle, txt, toc(t0));
end

vars = {'design', 'scenario', 'n', 'eff', 'sigTau', 'pSingle', 'estimator', ...
        'reject', 'reject_se', 'rejectStrat', 'nValid', 'bias', 'rmse', 'diffBias', 'diffRmse', 'coverage', ...
        'TnullMean', 'TnullP95', 'Tobs', 'pCal'};
R = cell2struct(rows, vars, 2);

save(fullfile(o.outDir, 'sim_results.mat'), 'R', 'o', 'G');
write_csv(R, vars, fullfile(o.outDir, 'sim_results.csv'));
write_tables(R, o, est);
if strcmp(o.design, 'calibrated'), write_calib_table(R, o, est); end
if o.bPlots
    try
        draw_figures(R, o, est);
    catch err
        warning('sim_framework:plots', 'Figures not drawn (%s). Results are saved.', err.message);
    end
end
fprintf('[sim] Done in %.1f min. Results in %s\n', toc(t0)/60, o.outDir);
end


%% =========================================================================
%  Generative model
%% =========================================================================
function G = default_gen()
G.T_RANGE   = [-10 5];   % analysis window (years from the event)
G.nq        = 150;       % query points on the common range
G.rangeP    = 10;        % coverage range: percentiles of subject medians
% Population trajectory of the reference subgroup (standardised scale):
% f0(t) = yHi - (yHi - yLo) * logistic((t - t50) / s)
G.yHi = 1.0;  G.yLo = -1.0;  G.t50 = 2;  G.s = 3;
% Positive subgroup, scenario 'diverge': the difference grows from 0 to
% -eff*SD around tDiv, with time scale sDiv.
G.tDiv = 0;  G.sDiv = 1.5;
% Between-subject (random intercept) and within-subject SDs.
G.sigB = 0.5;  G.sigE = 0.3;
% Location on the axis: subject centre ~ N(mu, sdC), clipped to cLim.
G.mu0 = -5;  G.mu1 = -1;  G.sdC = 3;  G.cLim = [-12 6];
% Visits: single-visit subjects with prob. pSingle; the others have
% 2 + Poisson(lambda * exp(gammaInf * b / sigB)) visits (max nMax),
% spaced 'gap' years with uniform jitter +-jit.
G.lambda = 2;  G.gammaInf = 0.5;  G.nMax = 12;  G.gap = 1;  G.jit = 0.25;
end

function y = f_true(t, g, C, G)
if ~isempty(C.cal)          % calibrated null: same global curve for both
    % Quadratic fitted inside the analysis window; held flat beyond a
    % 2-year margin so that distant visits do not follow an extrapolation.
    tc = t;
    if C.cal.robust, tc = min(max(t, G.T_RANGE(1) - 2), G.T_RANGE(2) + 2); end
    y = C.cal.c0 + C.cal.c1 * tc + C.cal.c2 * tc.^2;
    return
end
sd = sqrt(G.sigB^2 + G.sigE^2);
y = G.yHi - (G.yHi - G.yLo) ./ (1 + exp(-(t - G.t50) / G.s));
if g == 1
    switch C.scenario
        case 'level'
            y = y - C.eff * sd;
        case 'diverge'
            y = y - C.eff * sd ./ (1 + exp(-(t - G.tDiv) / G.sDiv));
    end
end
end

function D = generate_data(C, G)
% One replicate. D.id, D.g, D.t (observed axis), D.y; D.ttrue.
if isempty(C.cal)
    n = [C.n C.n];  mu = [G.mu0 G.mu1];  sdc = [G.sdC G.sdC];  cLim = G.cLim;
    ps = [C.pSingle C.pSingle];  lam = [G.lambda G.lambda];  gap = [G.gap G.gap];
    sigB = G.sigB;  sigE = G.sigE;
else
    K = C.cal;
    n = [K.n0 K.n1];  mu = [K.mu0 K.mu1];
    % Observed subject medians include the transition-age error; remove
    % the variance of its core component to get the true spread.
    sdc = sqrt(max([K.sd0 K.sd1].^2 - C.tau.core^2, 0.5^2));
    ps = [K.pSingle0 K.pSingle1];
    lam = max([K.nvMulti0 K.nvMulti1] - 2, 0.05);
    lam(isnan(lam)) = 0.05;
    gap = [K.gap0 K.gap1];  gap(isnan(gap)) = 1;
    sigB = K.sigB;  sigE = K.sigE;
    cLim = G.cLim;
    if K.robust, cLim = [-Inf Inf]; end   % real subjects also lie outside the window
end
nTot = sum(n);
g  = [zeros(n(1), 1); ones(n(2), 1)];
gi = g + 1;
cen = min(max(mu(gi)' + sdc(gi)' .* randn(nTot, 1), cLim(1)), cLim(2));
b   = sigB * randn(nTot, 1);
dlt = tau_error(C.tau, nTot);
if ~isempty(C.cal) && isfield(C.cal, 'noTau') && C.cal.noTau
    dlt(:) = 0;   % step 5 of the generator omitted; sdc keeps the core correction
end
single = rand(nTot, 1) < ps(gi)';
nv = ones(nTot, 1);
for i = find(~single)'
    nv(i) = min(2 + poisson_draw(lam(gi(i)) * exp(G.gammaInf * b(i) / max(sigB, eps))), G.nMax);
end
N = sum(nv);
id = zeros(N, 1); gg = id; tt = id; yy = id; to = id;
k = 0;
for i = 1:nTot
    j = (0:nv(i)-1)';
    ti = cen(i) + (j - (nv(i)-1)/2) * gap(gi(i)) + G.jit * (2*rand(nv(i), 1) - 1);
    yi = f_true(ti, g(i), C, G) + b(i) + sigE * randn(nv(i), 1);
    idx = k + (1:nv(i));
    id(idx) = i; gg(idx) = g(i); tt(idx) = ti; yy(idx) = yi; to(idx) = ti + dlt(i);
    k = k + nv(i);
end
D = struct('id', id, 'g', gg, 't', to, 'ttrue', tt, 'y', yy, ...
           'gSubj', g, 'nSubj', nTot, 'C', C);
end

function d = tau_error(T, n)
% Transition-age error: normal (core) or normal mixture (core + tail).
d = T.core * randn(n, 1);
if strcmp(T.model, 'mix')
    tail = rand(n, 1) < T.pTail;
    d(tail) = T.tail * randn(sum(tail), 1);
end
end

function T = default_tau_table()
% Fitted to tau_hat - clinical onset in converters (eval_dpm_onset_*):
% core = robust SD; tail chosen so that 0.85*core^2 + 0.15*tail^2 equals
% the classical SD^2 corrected for the visit interval.
T = struct('cohort',     {'ADNI',   'ADNI',    'NACC',   'NACC'}, ...
           'transition', {'CU2MCI', 'MCI2Dem', 'CU2MCI', 'MCI2Dem'}, ...
           'core',       {1.08,     1.11,      2.32,     1.92}, ...
           'tail',       {6.6,      5.1,       11.6,     6.9}, ...
           'pTail',      {0.15,     0.15,      0.15,     0.15});
end

function [cells, OBS] = calibrated_cells(o)
cells = struct('design', {}, 'scenario', {}, 'n', {}, 'eff', {}, ...
               'sigTau', {}, 'pSingle', {}, 'nReps', {}, 'cal', {}, 'tau', {});
files = cellstr(o.calibFiles);
if isempty(files), error('sim_framework:calib', 'calibFiles is required for the calibrated design.'); end
for f = 1:numel(files)
    Cs = read_csv_simple(files{f});
    for r = 1:numel(Cs.cohort)
        K = struct();
        fn = fieldnames(Cs);
        for k = 1:numel(fn)
            v = Cs.(fn{k});
            if iscell(v), K.(fn{k}) = v{r}; else, K.(fn{k}) = v(r); end
        end
        ti = find(strcmp({o.tauTable.cohort}, K.cohort) & strcmp({o.tauTable.transition}, K.transition), 1);
        if isempty(ti), error('sim_framework:tau', 'No tauTable entry for %s %s.', K.cohort, K.transition); end
        if any(isnan([o.tauTable(ti).core, o.tauTable(ti).tail]))
            error('sim_framework:tau', 'tauTable has NaN for %s %s. Fill in the transition-age error first.', K.cohort, K.transition);
        end
        tz = struct('model', 'mix', 'core', o.tauTable(ti).core, ...
                    'tail', o.tauTable(ti).tail, 'pTail', o.tauTable(ti).pTail);
        lab = sprintf('%s|%s|%s|%s', K.cohort, K.transition, K.analysis, K.var);
        K.robust = o.calRobust;
        K.noTau  = o.calNoTau;
        cells(end+1) = struct('design', 'calibrated', 'scenario', lab, 'n', K.n0 + K.n1, ...
            'eff', 0, 'sigTau', tz.core, 'pSingle', (K.pSingle0 + K.pSingle1) / 2, ...
            'nReps', o.nRepsCal, 'cal', K, 'tau', tz); %#ok<AGROW>
    end
end
OBS = struct('cohort', {}, 'transition', {}, 'analysis', {}, 'var', {}, 'bw', {}, 'T', {}, 'est', {});
for f = 1:numel(o.obsFiles)
    B = read_csv_simple(o.obsFiles{f});
    [~, nm] = fileparts(o.obsFiles{f});
    est = 'sw'; if numel(nm) > 3 && strcmp(nm(end-2:end), '_ll'), est = 'swll'; end
    for r = 1:numel(B.cohort)
        OBS(end+1) = struct('cohort', B.cohort{r}, 'transition', B.transition{r}, ...
            'analysis', B.analysis{r}, 'var', B.var{r}, 'bw', B.bw(r), 'T', B.T_obs(r), 'est', est); %#ok<AGROW>
    end
end
end

function T = lookup_obs(OBS, C, e, bw)
T = NaN;
if isempty(OBS) || ~strcmp(C.design, 'calibrated'), return; end
if strncmp(e, 'swll_bw', 7), est = 'swll'; elseif strncmp(e, 'sw_bw', 5), est = 'sw'; else, return; end
K = C.cal;
m = find(strcmp({OBS.cohort}, K.cohort) & strcmp({OBS.transition}, K.transition) & ...
         strcmp({OBS.analysis}, K.analysis) & strcmp({OBS.var}, K.var) & ...
         abs([OBS.bw] - bw) < 1e-9 & strcmp({OBS.est}, est), 1);
if ~isempty(m), T = OBS(m).T; end
end

function S = read_csv_simple(file)
% Minimal CSV reader (header + rows; numeric columns converted).
fid = fopen(file, 'r');
if fid < 0, error('sim_framework:csv', 'Cannot open %s', file); end
hdr = strtrim(strsplit(fgetl(fid), ','));
hdr = strrep(hdr, '"', '');
rows = {};
while true
    l = fgetl(fid);
    if ~ischar(l), break; end
    if isempty(strtrim(l)), continue; end
    rows(end+1, :) = strrep(strtrim(strsplit(l, ',')), '"', ''); %#ok<AGROW>
end
fclose(fid);
S = struct();
for k = 1:numel(hdr)
    col = rows(:, k);
    num = str2double(col);
    nm = matlab.lang.makeValidName(hdr{k});
    if all(~isnan(num) | strcmpi(col, 'NaN') | cellfun(@isempty, col))
        S.(nm) = num;
    else
        S.(nm) = col;
    end
end
end

function k = poisson_draw(lam)
% Knuth's algorithm (lam is small).
L = exp(-lam); k = 0; q = rand;
while q > L
    k = k + 1; q = q * rand;
end
end


%% =========================================================================
%  Analysis of one replicate
%% =========================================================================
function est = estimator_names(o)
est = arrayfun(@(b) sprintf('sw_bw%g', b), o.bwGrid, 'UniformOutput', false);
est = [est, arrayfun(@(b) sprintf('swll_bw%g', b), o.bwGrid, 'UniformOutput', false), {'nw', 'poly'}];
end

function M = analyse(D, G, o, est, doBoot)
ne = numel(est);
M.p = nan(2, ne); M.bias = nan(1, ne); M.rmse = M.bias; M.dbias = M.bias; M.drmse = M.bias;
M.cover = M.bias; M.T = M.bias;

in0 = D.gSubj == 0;
tr0 = coverage_range(D, find(in0), G);
tr1 = coverage_range(D, find(~in0), G);
lo = max(tr0(1), tr1(1)); hi = min(tr0(2), tr1(2));
if hi <= lo, return; end
xq = linspace(lo, hi, G.nq);
f0 = f_true(xq, 0, D.C, G);  f1 = f_true(xq, 1, D.C, G);

% Same permutations for every estimator. Two schemes:
%   free  : labels shuffled among all subjects (manuscript, Sec. 3.3.2);
%   strat : labels shuffled within strata of the subject-median disease
%           time, so that each permuted subgroup keeps the location of the
%           original one on the axis (conditional permutation test).
n0 = sum(in0); nS = D.nSubj;
P0 = false(nS, o.nPerm);
for b = 1:o.nPerm
    idx = randperm(nS); P0(idx(1:n0), b) = true;
end
P1 = ~P0;
Q0 = strat_perms(D, in0, o.nPerm, o.nStrata);
Q1 = ~Q0;

for k = 1:ne
    e = est{k};
    if strncmp(e, 'sw_bw', 5) || strcmp(e, 'nw')
        if strcmp(e, 'nw'), bw = o.bw; mode = 'visit';
        else, bw = sscanf(e, 'sw_bw%f'); mode = 'subject'; end
        [A, B] = kernel_terms(D, xq, bw, mode);
        cf = @(Msk) curve(A, B, Msk);
    elseif strncmp(e, 'swll_bw', 7)
        L = ll_terms(D, xq, sscanf(e, 'swll_bw%f'));
        cf = @(Msk) ll_curve(L, Msk);
    else
        S = poly_terms(D, G);
        cf = @(Msk) poly_curve(S, Msk, xq);
    end
    y0 = cf(in0);  y1 = cf(~in0);
    [M.p(1, k), ok, M.T(k)] = perm_p(y0, y1, cf(P0), cf(P1));
    if ~ok, continue; end
    M.p(2, k) = perm_p(y0, y1, cf(Q0), cf(Q1));
    M.bias(k)  = mean([y0 - f0, y1 - f1], 'omitnan');
    M.rmse(k)  = sqrt(mean([y0 - f0, y1 - f1].^2, 'omitnan'));
    d = (y1 - y0) - (f1 - f0);
    M.dbias(k) = mean(d, 'omitnan');
    M.drmse(k) = sqrt(mean(d.^2, 'omitnan'));
    % Bootstrap bands for the main sw-LOESS bandwidth.
    if doBoot && (strcmp(e, sprintf('sw_bw%g', o.bw)) || strcmp(e, sprintf('swll_bw%g', o.bw)))
        M.cover(k) = mean([boot_cover(cf, D.nSubj, find(in0), f0, o.nBoot), ...
                           boot_cover(cf, D.nSubj, find(~in0), f1, o.nBoot)], 'omitnan');
    end
end
end

function tr = coverage_range(D, subj, G)
% Percentiles of the subject-median disease time, clipped to the window;
% full window if the range covers >= 90% of it (as in plot_tau_AV.m).
med = zeros(numel(subj), 1);
for k = 1:numel(subj)
    med(k) = median(D.t(D.id == subj(k)));
end
tr = G.T_RANGE;
if numel(med) < 3, return; end
lo = max(prctile(med, G.rangeP), G.T_RANGE(1));
hi = min(prctile(med, 100 - G.rangeP), G.T_RANGE(2));
if hi <= lo, return; end
if (hi - lo) < 0.90 * diff(G.T_RANGE), tr = [lo hi]; end
end

function [A, B] = kernel_terms(D, xq, bw, mode)
% Per-subject terms on the query grid, so that a subgroup curve is
% sum(B(subset,:)) ./ sum(A(subset,:)).
%   'subject' (sw-LOESS): A = w_i(t) = max_j k_ij(t), B = w_i(t) * ybar_i(t)
%   'visit'   (classical NW): A = sum_j k_ij(t),      B = sum_j k_ij(t) y_ij
nS = D.nSubj; nq = numel(xq);
A = zeros(nS, nq); B = A;
for i = 1:nS
    m = D.id == i;
    K = exp(-0.5 * ((xq - D.t(m)) / bw).^2);     % visits x nq
    sk = sum(K, 1);  sky = sum(K .* D.y(m), 1);
    if strcmp(mode, 'subject')
        w = max(K, [], 1);
        A(i, :) = w;  B(i, :) = w .* sky ./ max(sk, realmin);
    else
        A(i, :) = sk; B(i, :) = sky;
    end
end
end

function Y = curve(A, B, Msk)
% Msk: nS x k logical (one column per subset). Returns k x nq.
% Fewer than 3 subjects, or no weight, gives NaN (as eval_curve).
Mf = double(Msk);
num = Mf' * B;  den = Mf' * A;
Y = num ./ den;
Y(den <= 1e-12) = NaN;
Y(sum(Mf > 0, 1)' < 3, :) = NaN;
end

function Q0 = strat_perms(D, in0, nPerm, nStrata)
% Permutations within strata of the subject-median disease time. Within
% each stratum the labels are shuffled, so the number of reference
% subjects per stratum is preserved.
nS = D.nSubj;
med = zeros(nS, 1);
for i = 1:nS, med(i) = median(D.t(D.id == i)); end
edges = prctile(med, linspace(0, 100, nStrata + 1));
edges(1) = -Inf; edges(end) = Inf;
st = zeros(nS, 1);
for k = 1:nStrata, st(med > edges(k) & med <= edges(k+1)) = k; end
Q0 = false(nS, nPerm);
for k = 1:nStrata
    m = find(st == k); lab = in0(m);
    for b = 1:nPerm
        Q0(m, b) = lab(randperm(numel(m)));
    end
end
end

function L = ll_terms(D, xq, bw)
% Subject-weighted local-linear smoother. At each t, a weighted line is
% fitted with visit weights a_ij = w_i(t) k_ij(t) / sum_j k_ij(t), so each
% subject keeps the total weight w_i(t) = max_j k_ij(t) of sw-LOESS.
% Per-subject moments are stored so subset fits use matrix products.
nS = D.nSubj; nq = numel(xq);
L.s0 = zeros(nS, nq); L.s1 = L.s0; L.s2 = L.s0; L.t0 = L.s0; L.t1 = L.s0;
for i = 1:nS
    m = D.id == i;
    dt = D.t(m) - xq;                              % visits x nq
    K = exp(-0.5 * (dt / bw).^2);
    a = K .* (max(K, [], 1) ./ max(sum(K, 1), realmin));
    y = D.y(m);
    L.s0(i, :) = sum(a, 1);          L.s1(i, :) = sum(a .* dt, 1);
    L.s2(i, :) = sum(a .* dt.^2, 1); L.t0(i, :) = sum(a .* y, 1);
    L.t1(i, :) = sum(a .* dt .* y, 1);
end
end

function Y = ll_curve(L, Msk)
Mf = double(Msk)';
S0 = Mf * L.s0; S1 = Mf * L.s1; S2 = Mf * L.s2; T0 = Mf * L.t0; T1 = Mf * L.t1;
den = S0 .* S2 - S1.^2;
Y = (S2 .* T0 - S1 .* T1) ./ den;
Y(den <= 1e-10 * max(S0.^2, realmin)) = NaN;
Y(sum(Mf > 0, 2) < 3, :) = NaN;
end

function S = poly_terms(D, G)
% Per-subject sufficient statistics of a quadratic fit on the visits that
% fall inside the analysis window. Time is scaled by 1/5 for conditioning.
nS = D.nSubj;
S.m = zeros(nS, 5); S.v = zeros(nS, 3);
in = D.t >= G.T_RANGE(1) & D.t <= G.T_RANGE(2);
u = D.t / 5;
for i = 1:nS
    m = D.id == i & in;
    ui = u(m); yi = D.y(m);
    S.m(i, :) = [numel(ui), sum(ui), sum(ui.^2), sum(ui.^3), sum(ui.^4)];
    S.v(i, :) = [sum(yi), sum(ui .* yi), sum(ui.^2 .* yi)];
end
end

function Y = poly_curve(S, Msk, xq)
Mf = double(Msk);
Mm = Mf' * S.m;  Vv = Mf' * S.v;
u = xq / 5;  X = [ones(size(u)); u; u.^2];
Y = nan(size(Msk, 2), numel(xq));
for k = 1:size(Msk, 2)
    m = Mm(k, :);
    if m(1) < 4, continue; end
    XtX = [m(1) m(2) m(3); m(2) m(3) m(4); m(3) m(4) m(5)];
    if rcond(XtX) < 1e-12, continue; end
    beta = XtX \ Vv(k, :)';
    Y(k, :) = beta' * X;
end
end

function [p, ok, Tobs] = perm_p(y0, y1, Y0p, Y1p)
okq = ~isnan(y0) & ~isnan(y1);
ok = sum(okq) >= 5;
Tobs = NaN;
if ~ok, p = NaN; return; end
Tobs = mean(abs(y1(okq) - y0(okq)));
Dp = abs(Y1p - Y0p);
valid = ~isnan(Dp);
Dp(~valid) = 0;
Tp = sum(Dp, 2) ./ max(sum(valid, 2), 1);
Tp(sum(valid, 2) < 3) = 0;        % 'zero' rule of the original scripts
p = (1 + sum(Tp >= Tobs)) / (1 + numel(Tp));
end

function cv = boot_cover(cf, nS, subj, ftrue, nBoot)
% Pointwise 95% percentile band from a subject bootstrap within the
% subgroup; returns the fraction of query points covered by the band.
% Resampled subjects enter as integer weights (counts).
n = numel(subj);
Cn = zeros(nS, nBoot);
for b = 1:nBoot
    Cn(:, b) = accumarray(subj(randi(n, n, 1)), 1, [nS 1]);
end
Yb = cf(Cn);
lo = prctile(Yb, 2.5, 1);  hi = prctile(Yb, 97.5, 1);
ok = ~isnan(lo) & ~isnan(hi);
cv = mean(ftrue(ok) >= lo(ok) & ftrue(ok) <= hi(ok));
end


%% =========================================================================
%  Output
%% =========================================================================
function v = ifelse(c, a, b)
if c, v = a; else, v = b; end
end

function write_csv(R, vars, file)
fid = fopen(file, 'w');
fprintf(fid, '%s\n', strjoin(vars, ','));
for r = 1:numel(R)
    x = R(r);
    fprintf(fid, '%s,%s,%d,%.3f,%.2f,%.2f,%s,%.4f,%.4f,%.4f,%d,%.5f,%.5f,%.5f,%.5f,%.4f,%.5f,%.5f,%.5f,%.4f\n', ...
        x.design, x.scenario, x.n, x.eff, x.sigTau, x.pSingle, x.estimator, ...
        x.reject, x.reject_se, x.rejectStrat, x.nValid, x.bias, x.rmse, x.diffBias, x.diffRmse, x.coverage, ...
        x.TnullMean, x.TnullP95, x.Tobs, x.pCal);
end
fclose(fid);
end

function s = est_label(e)
if strncmp(e, 'sw_bw', 5), s = sprintf('local-constant \\textit{sw}-LOESS, $bw=%s$', e(6:end));
elseif strncmp(e, 'swll_bw', 7), s = sprintf('local-linear \\textit{sw}-LOESS, $bw=%s$', e(8:end));
elseif strcmp(e, 'nw'), s = 'visit-weighted NW';
else, s = 'global polynomial'; end
end

function s = dec(x, n)
if isnan(x), s = '--'; return; end
s = sprintf(['%.' num2str(n) 'f'], x);
end

function write_tables(R, o, est)
% LaTeX tables (Supplementary Material style).
d = {R.design}; sc = {R.scenario}; es = {R.estimator};
n = [R.n]; ef = [R.eff]; st = [R.sigTau]; ps = [R.pSingle];

% --- type I
if any(strcmp(d, 'typeI'))
    cfg = unique([st(strcmp(d,'typeI'))' ps(strcmp(d,'typeI'))'], 'rows');
    ns = unique(n(strcmp(d, 'typeI')));
    fid = fopen(fullfile(o.outDir, 'sim_typeI.tex'), 'w');
    fprintf(fid, '%% Auto-generated by sim_framework.m\n\\begin{table}[htbp]\n\\centering\n');
    fprintf(fid, ['\\caption{Type~I error of the permutation test (proportion of replicates with $p<0.05$) ' ...
        'when the two subgroups share the same trajectory. $n$: subjects per subgroup. ' ...
        '$\\sigma_\\tau$: SD of the transition-age error (years). ' ...
        '$p_1$: proportion of subjects with a single visit. %d replicates per cell and %d permutations. ' ...
        'The nominal value is 0.05, with a Monte Carlo error of about %s.}\n'], ...
        o.nRepsNull, o.nPerm, dec(sqrt(0.05*0.95/o.nRepsNull), 3));
    fprintf(fid, '\\label{tab:S_sim_typeI}\n\\scriptsize\n\\setlength{\\tabcolsep}{4pt}\n');
    fprintf(fid, '\\begin{tabular}{@{}lc%s@{}}\n\\toprule\n', repmat('c', 1, size(cfg,1)*numel(ns)));
    fprintf(fid, ' & ');
    for c = 1:size(cfg, 1)
        fprintf(fid, ' & \\multicolumn{%d}{c}{$\\sigma_\\tau=%s$, $p_1=%s$}', numel(ns), dec(cfg(c,1),1), dec(cfg(c,2),2));
    end
    fprintf(fid, ' \\\\\n\\textbf{Estimator} & ');
    for c = 1:size(cfg, 1), fprintf(fid, '%s', sprintf(' & $n=%d$', ns)); end
    fprintf(fid, ' \\\\\n\\midrule\n');
    for k = 1:numel(est)
        fprintf(fid, '%s & ', est_label(est{k}));
        for c = 1:size(cfg, 1)
            for j = ns
                r = find(strcmp(d,'typeI') & strcmp(es, est{k}) & n == j & st == cfg(c,1) & ps == cfg(c,2), 1);
                fprintf(fid, ' & %s', dec(R(r).reject, 3));
            end
        end
        fprintf(fid, ' \\\\\n');
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
    fclose(fid);
end

% --- power and bias
if any(strcmp(d, 'power'))
    scs = unique(sc(strcmp(d, 'power')));
    ns = unique(n(strcmp(d, 'power'))); es_ = unique(ef(strcmp(d, 'power')));
    names = {'sim_power.tex', 'sim_bias.tex'};
    for t = 1:2
        fid = fopen(fullfile(o.outDir, names{t}), 'w');
        fprintf(fid, '%% Auto-generated by sim_framework.m\n\\begin{table}[htbp]\n\\centering\n');
        if t == 1
            fprintf(fid, ['\\caption{Power of the permutation test (proportion of replicates with $p<0.05$). ' ...
                'Level: the positive subgroup is shifted along the whole axis. Divergence: the difference grows ' ...
                'around the event. Effect size in SD of the marker. ' ...
                '$\\sigma_\\tau=%s$~years, $p_1=%s$, %d replicates per cell.}\n'], dec(o.sigTau,1), dec(o.pSingle,2), o.nRepsAlt);
            fprintf(fid, '\\label{tab:S_sim_power}\n');
        else
            fprintf(fid, ['\\caption{Mean bias of the difference between curves, $(\\hat{y}_1-\\hat{y}_0)-(f_1-f_0)$, ' ...
                'on the common range. A positive value means that the estimated difference is smaller than the true one ' ...
                '(attenuation). Same replicates as Table~\\ref{tab:S_sim_power}.}\n']);
            fprintf(fid, '\\label{tab:S_sim_bias}\n');
        end
        fprintf(fid, '\\scriptsize\n\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{@{}ll%s@{}}\n\\toprule\n', ...
            repmat('c', 1, numel(es_) * numel(scs)));
        fprintf(fid, ' & ');
        for s = 1:numel(scs)
            lab = ifelse(strcmp(scs{s}, 'level'), 'Level', 'Divergence');
            fprintf(fid, ' & \\multicolumn{%d}{c}{%s}', numel(es_), lab);
        end
        fprintf(fid, ' \\\\\n$n$ & \\textbf{Estimator}');
        for s = 1:numel(scs), fprintf(fid, '%s', sprintf(' & %s', strjoin(arrayfun(@(x) dec(x,2), es_, 'UniformOutput', false), ' & '))); end
        fprintf(fid, ' \\\\\n\\midrule\n');
        for j = ns
            for k = 1:numel(est)
                fprintf(fid, '%s & %s', ifelse(k == 1, sprintf('%d', j), ''), est_label(est{k}));
                for s = 1:numel(scs)
                    for e = es_
                        r = find(strcmp(d,'power') & strcmp(sc, scs{s}) & strcmp(es, est{k}) & n == j & ef == e, 1);
                        if t == 1, v = dec(R(r).reject, 2); else, v = dec(R(r).diffBias, 3); end
                        fprintf(fid, ' & %s', v);
                    end
                end
                fprintf(fid, ' \\\\\n');
            end
            if j ~= ns(end), fprintf(fid, '\\midrule\n'); end
        end
        fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
        fclose(fid);
    end

    % --- coverage
    fid = fopen(fullfile(o.outDir, 'sim_coverage.tex'), 'w');
    fprintf(fid, '%% Auto-generated by sim_framework.m\n\\begin{table}[htbp]\n\\centering\n');
    fprintf(fid, ['\\caption{Pointwise coverage of the 95\\,\\%% \\textit{bootstrap} bands of local-constant \\textit{sw}-LOESS ' ...
        '($bw=%s$): proportion of the common range where the band contains the true curve. ' ...
        '%d resamples per replicate.}\n'], dec(o.bw, 0), o.nBoot);
    fprintf(fid, '\\label{tab:S_sim_cov}\n\\scriptsize\n\\begin{tabular}{@{}l%s@{}}\n\\toprule\n', ...
        repmat('c', 1, numel(es_) * numel(scs)));
    for s = 1:numel(scs)
        lab = ifelse(strcmp(scs{s}, 'level'), 'Level', 'Divergence');
        fprintf(fid, ' & \\multicolumn{%d}{c}{%s}', numel(es_), lab);
    end
    fprintf(fid, ' \\\\\n$n$');
    for s = 1:numel(scs), fprintf(fid, '%s', sprintf(' & %s', strjoin(arrayfun(@(x) dec(x,2), es_, 'UniformOutput', false), ' & '))); end
    fprintf(fid, ' \\\\\n\\midrule\n');
    em = sprintf('sw_bw%g', o.bw);
    for j = ns
        fprintf(fid, '%d', j);
        for s = 1:numel(scs)
            for e = es_
                r = find(strcmp(d,'power') & strcmp(sc, scs{s}) & strcmp(es, em) & n == j & ef == e, 1);
                fprintf(fid, ' & %s', dec(R(r).coverage, 3));
            end
        end
        fprintf(fid, ' \\\\\n');
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
    fclose(fid);
end
end

function write_calib_table(R, o, est)
% One row per real comparison and marker: type I error with the realistic
% design and transition-age error, mean and 95th percentile of the null
% statistic, observed statistic and calibrated p-value.
idx = find(strcmp({R.design}, 'calibrated'));
if isempty(idx), return; end
labs = unique({R(idx).scenario}, 'stable');
fid = fopen(fullfile(o.outDir, 'sim_calibrated.tex'), 'w');
fprintf(fid, '%% Auto-generated by sim_framework.m\n\\begin{table}[htbp]\n\\centering\n');
fprintf(fid, ['\\caption{Calibrated simulation of each real comparison. Null hypothesis with the ' ...
    'same curve in the two subgroups, the observed design and the transition-age error ' ...
    'of each cohort and transition. Type~I: proportion of replicates with $p<0.05$. ' ...
    '$\\bar{T}_0$ and $T_{0,95}$: mean and 95th percentile of the statistic under the null. ' ...
    '$T_{\\rm obs}$: observed statistic. $p_{\\rm cal}$: calibrated $p$-value, $(1+\\#\\{T_0\\geq T_{\\rm obs}\\})/(1+R)$. ' ...
    '%d replicates per comparison, $bw=%s$.}\n'], o.nRepsCal, dec(o.bw, 0));
fprintf(fid, '\\label{tab:S_sim_cal}\n\\scriptsize\n\\setlength{\\tabcolsep}{3pt}\n');
fprintf(fid, '\\begin{tabular}{@{}llll ccccc ccccc c@{}}\n\\toprule\n');
fprintf(fid, [' & & & & \\multicolumn{5}{c}{local-linear \\textit{sw}-LOESS} & \\multicolumn{5}{c}{local-constant \\textit{sw}-LOESS} & Polynomial \\\\\n' ...
    '\\cmidrule(lr){5-9}\\cmidrule(lr){10-14}\\cmidrule(lr){15-15}\n' ...
    'Cohort & Trans. & Comp. & Marker & Type~I & $\\bar{T}_0$ & $T_{0,95}$ & $T_{\\rm obs}$ & $p_{\\rm cal}$ ' ...
    '& Type~I & $\\bar{T}_0$ & $T_{0,95}$ & $T_{\\rm obs}$ & $p_{\\rm cal}$ & Type~I \\\\\n\\midrule\n']);
el = sprintf('swll_bw%g', o.bw);  ec = sprintf('sw_bw%g', o.bw);
for k = 1:numel(labs)
    parts = strsplit(labs{k}, '|');
    rl = idx(strcmp({R(idx).scenario}, labs{k}) & strcmp({R(idx).estimator}, el));
    rc = idx(strcmp({R(idx).scenario}, labs{k}) & strcmp({R(idx).estimator}, ec));
    rp = idx(strcmp({R(idx).scenario}, labs{k}) & strcmp({R(idx).estimator}, 'poly'));
    fprintf(fid, '%s & %s & %s & %s', parts{1}, strrep(parts{2}, '2', '$\to$'), ...
        strrep(strrep(parts{3}, '_', ' '), '+', '$+$'), parts{4});
    for rr = [rl rc]
        x = R(rr);
        fprintf(fid, ' & %s & %s & %s & %s & %s', dec(x.reject, 2), dec(x.TnullMean, 3), ...
            dec(x.TnullP95, 3), dec(x.Tobs, 3), dec(x.pCal, 3));
    end
    fprintf(fid, ' & %s \\\\\n', dec(R(rp).reject, 2));
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
fclose(fid);
end

function draw_figures(R, o, est)
d = {R.design}; sc = {R.scenario}; es = {R.estimator};
n = [R.n]; ef = [R.eff];
if any(strcmp(d, 'power'))
    scs = unique(sc(strcmp(d, 'power'))); ns = unique(n(strcmp(d, 'power')));
    es_ = unique(ef(strcmp(d, 'power')));
    fh = figure('Color', 'w', 'Position', [100 100 300*numel(ns) 260*numel(scs)]);
    for s = 1:numel(scs)
        for j = 1:numel(ns)
            subplot(numel(scs), numel(ns), (s-1)*numel(ns) + j); hold on; box on; grid on;
            for k = 1:numel(est)
                y = arrayfun(@(e) R(find(strcmp(d,'power') & strcmp(sc,scs{s}) & strcmp(es,est{k}) & n==ns(j) & ef==e, 1)).reject, es_);
                plot(es_, y, '-o', 'LineWidth', 1.2, 'DisplayName', strrep(est{k}, '_', ' '));
            end
            ylim([0 1]); xlabel('Effect (SD)'); ylabel('Power');
            title(sprintf('%s, n = %d', scs{s}, ns(j)));
            if s == 1 && j == numel(ns), legend('Location', 'southeast'); end
        end
    end
    print(fh, fullfile(o.outDir, 'sim_power.png'), '-dpng', '-r150');
end
if any(strcmp(d, 'typeI'))
    idx = find(strcmp(d, 'typeI'));
    lab = arrayfun(@(r) sprintf('n=%d st=%.1f p1=%.2f', R(r).n, R(r).sigTau, R(r).pSingle), idx, 'UniformOutput', false);
    ul = unique(lab, 'stable');
    [~, gi] = ismember(lab, ul);  gi = gi(:);
    eidx = es(idx);
    Y = nan(numel(ul), numel(est));
    for k = 1:numel(est)
        for u = 1:numel(ul)
            m = gi == u & reshape(strcmp(eidx, est{k}), [], 1);
            r = idx(m);
            if ~isempty(r), Y(u, k) = R(r(1)).reject; end
        end
    end
    fh = figure('Color', 'w', 'Position', [100 100 900 350]);
    bar(Y); hold on; plot([0.5 numel(ul)+0.5], [0.05 0.05], 'k--');
    set(gca, 'XTick', 1:numel(ul), 'XTickLabel', ul, 'XTickLabelRotation', 30);
    ylabel('Type I error'); legend(strrep(est, '_', ' '), 'Location', 'northwest'); box on;
    print(fh, fullfile(o.outDir, 'sim_typeI.png'), '-dpng', '-r150');
end
end
