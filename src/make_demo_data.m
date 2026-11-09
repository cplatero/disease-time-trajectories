function make_demo_data(dataDir, tauDir, seed)
%MAKE_DEMO_DATA  Synthetic cohorts with the layout of the real CSV files.
%
%   make_demo_data('data/demo', 'data/demo/tau')
%
%   Writes ADNI_<TR>.csv, NACC_<TR>.csv (visits) and tau_<COHORT>_<TR>.csv
%   (10 transition-age estimates per subject) for TR = CU2MCI, MCI2Dem.
%   They let anyone run every experiment and the tutorial without access
%   to ADNI or NACC. The numbers have NO clinical meaning: only the
%   structure of the data is realistic.
%
%   Generative model (per subject i):
%     clinical group: stable (1), converter (2) or next stage (3), with the
%       true transition age tau_i after, within or before follow-up;
%     amyloid status A+ (more frequent in groups 2-3) and baseline WMH
%       (lognormal; V+ if > 10 ml);
%     visits every year, 1 + Poisson(3) visits, age at baseline ~ N(72, 6);
%     MEM, EF: sigmoid decline on the true disease time t = age - tau_i,
%       lower in A+ (memory) and in V+ (executive function), plus a
%       subject intercept and visit noise;
%     Centiloids and WMH: one or two measures per subject;
%     transition-age estimates: tau_i + delta_i + model noise, with delta_i
%       from the 85/15 normal mixture of the paper.
%
%   Authorship: Carlos Platero (UPM), with the collaboration of Claude
%   (Anthropic).

if nargin < 1 || isempty(dataDir), dataDir = fullfile('data', 'demo'); end
if nargin < 2 || isempty(tauDir),  tauDir  = fullfile(dataDir, 'tau'); end
if nargin < 3 || isempty(seed),    seed    = 1; end
if ~exist(dataDir, 'dir'), mkdir(dataDir); end
if ~exist(tauDir, 'dir'),  mkdir(tauDir);  end
rng(seed, 'twister');

COH = {'ADNI', 'NACC'};
NS  = [500, 650];                  % subjects per transition
TAU = [1.1 6.0; 2.0 9.0];          % core and tail SD of delta_i per cohort
TR  = {'CU2MCI', 'MCI2Dem'};
for c = 1:2
    for t = 1:2
        [V, S] = simulate(NS(c), TAU(c, :), t, c * 10 + t);
        write_csv(fullfile(dataDir, sprintf('%s_%s.csv', COH{c}, TR{t})), V, ...
            {'ID', 'AGE_visit', 'ClinGroup', 'MEM', 'EF', 'CENTILOIDS', 'WMH', ...
             'A_bl', 'CENTILOIDS_bl', 'WMH_bl', 'onset_age_clin'});
        write_csv(fullfile(tauDir, sprintf('tau_%s_%s.csv', COH{c}, TR{t})), S, ...
            [{'ID'}, arrayfun(@(j) sprintf('tau_%d', j), 1:10, 'UniformOutput', false)]);
        fprintf('[demo] %s %s: %d subjects, %d visits\n', COH{c}, TR{t}, size(S, 1), size(V, 1));
    end
end
end


function [V, S] = simulate(N, tauSD, tr, idBase)
u = rand(N, 1);
grp = 1 + (u >= 0.50) + (u >= 0.75);                    % 1: 50 %, 2 and 3: 25 % each
pA  = [0.30 0.60 0.70];
A   = rand(N, 1) < pA(grp)';
WMHbl = exp(log(5) + 0.9 * randn(N, 1) + 0.3 * (grp > 1));
Vp  = WMHbl > 10;
age0 = 72 + 6 * randn(N, 1);
nv  = 1 + poissrnd_(3, N);
fu  = nv - 1;                                           % years of follow-up
tau = zeros(N, 1);
for i = 1:N
    switch grp(i)
        case 1, tau(i) = age0(i) + fu(i) + 1 + 8 * rand;          % after follow-up
        case 2, tau(i) = age0(i) + 0.5 + max(fu(i) - 0.5, 0.5) * rand;  % during follow-up
        case 3, tau(i) = age0(i) - 0.5 - 5 * rand;                % before baseline
    end
end
% converters need at least two visits
nv(grp == 2) = max(nv(grp == 2), 2);  fu = nv - 1;
shift = 2 * (tr == 2);                    % MCI -> dementia: later stage, faster decline
bS = 0.30 * randn(N, 1);  bE = 0.30 * randn(N, 1);
CLbl = A .* (45 + 25 * randn(N, 1)) + ~A .* (2 + 6 * randn(N, 1));
CLbl = max(CLbl, -15);

rows = zeros(sum(nv), 11);  r = 0;
for i = 1:N
    nPET = 1 + (rand < 0.3);  nMRI = 1 + (rand < 0.4);
    onset = NaN;
    for j = 0:nv(i) - 1
        age = age0(i) + j;
        tt  = age - tau(i);
        mem = 0.9 - 1.6 ./ (1 + exp(-(tt + 1 - shift) / 2.0)) - 0.25 * A(i) .* (tt > -6) ...
              + bS(i) + 0.20 * randn;
        ef  = 0.8 - 1.3 ./ (1 + exp(-(tt + 0.5 - shift) / 2.5)) - 0.20 * Vp(i) ...
              + bE(i) + 0.22 * randn;
        cl = NaN;  wm = NaN;
        if j < nPET * 2 && mod(j, 2) == 0, cl = CLbl(i) + 3 * A(i) * j + 4 * randn; end
        if j < nMRI * 2 && mod(j, 2) == 0, wm = WMHbl(i) * exp(0.05 * j + 0.08 * randn); end
        if grp(i) == 2 && isnan(onset) && age >= tau(i), onset = age; end
        r = r + 1;
        rows(r, :) = [idBase * 10000 + i, age, grp(i), mem, ef, cl, wm, double(A(i)), ...
                      CLbl(i), WMHbl(i), NaN];
    end
    rows(r - nv(i) + 1:r, 11) = onset;
end
% missing baseline biomarkers
noPET = rand(N, 1) < 0.25;  noMRI = rand(N, 1) < 0.15;
ids = rows(:, 1) - idBase * 10000;
rows(noPET(ids), [6 8 9]) = NaN;
rows(noMRI(ids), [7 10])  = NaN;
V = rows;

% transition-age estimates: common error delta_i plus model noise
core = tauSD(1);  tail = tauSD(2);
delta = core * randn(N, 1);
isT = rand(N, 1) < 0.15;
delta(isT) = tail * randn(sum(isT), 1);
E = tau + delta + 0.3 * randn(N, 10);
S = [idBase * 10000 + (1:N)', E];
end


function k = poissrnd_(lam, n)
% Poisson draws without the Statistics Toolbox (Knuth).
k = zeros(n, 1);
L = exp(-lam);
for i = 1:n
    p = 1;  x = -1;
    while p > L
        x = x + 1;  p = p * rand;
    end
    k(i) = x;
end
end


function write_csv(file, M, names)
fid = fopen(file, 'w');
fprintf(fid, '%s\n', strjoin(names, ','));
fmt = [strjoin(repmat({'%.6g'}, 1, numel(names)), ','), '\n'];
fprintf(fid, fmt, M');
fclose(fid);
end
