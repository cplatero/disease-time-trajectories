function C = config(mode)
%CONFIG  Paths and run settings of every experiment.
%
%   C = config('real')   cohort data in data/, results in results/real/
%   C = config('demo')   synthetic data in data/demo/, results in results/demo/
%
%   The methodological constants of the paper (disease-time window
%   [-10, 5] years, coverage percentiles 10-90, 150 query points, subject
%   bootstrap of 300 resamples, 1000 permutations of the kinematics)
%   are fixed inside bw_sensitivity.m and plot_kinematics.m, so that the
%   published numbers cannot be changed by accident.
%
%   Authorship: Carlos Platero (UPM), with the collaboration of Claude
%   (Anthropic).

if nargin < 1 || isempty(mode), mode = 'real'; end
mode = lower(char(mode));
if ~any(strcmp(mode, {'real', 'demo'}))
    error('config:mode', 'mode must be ''real'' or ''demo''.');
end

C.mode = mode;
C.root = fileparts(mfilename('fullpath'));
if strcmp(mode, 'real')
    C.dataDir = fullfile(C.root, 'data');
    C.tauDir  = fullfile(C.root, 'data', 'tau');
else
    C.dataDir = fullfile(C.root, 'data', 'demo');
    C.tauDir  = fullfile(C.root, 'data', 'demo', 'tau');
end
C.resultsDir = fullfile(C.root, 'results', mode);
C.cohorts    = {'ADNI', 'NACC'};

% Stratification thresholds (Sec. 2 of the paper)
C.CL_THRESHOLD  = 20;      % A+ : Centiloids_bl > 20 (applied when the CSV files are written)
C.WMH_THRESHOLD = 10;      % V+ : WMH_bl > 10 ml

% Comparisons. Table 4 uses A+ vs A- and V+ vs V- (24 comparisons).
% Fig. 3 and Table S9 also include V+ vs V- within A+ and within A-
% (48 comparisons in total). The A x V comparisons share the random
% stream of the calibration, so with {'A', 'V', 'AxV_A+', 'AxV_A-'} every
% number of the paper is reproduced exactly. With {'A', 'V'} the run is
% about twice as fast, Fig. 3 shows 24 comparisons, and the p_cal of
% MCI -> dementia change only by Monte Carlo error.
C.ANALYSES = {'A', 'V', 'AxV_A+', 'AxV_A-'};

% Bandwidth of sw-LOESS and grid of the sensitivity analysis (years)
C.BW      = 3;
C.BW_GRID = [2 3 4 6 8];
C.SEED    = 2026;

% Transition-age error (Table 1). 'paper' uses the published values, so
% the calibrated p-values are reproduced exactly; 'computed' uses the
% values estimated by exp2_transition_age_error.m from the CSV files.
C.TAU_SOURCE = 'paper';
C.TAU_PAPER  = struct( ...
    'cohort',     {'ADNI',   'ADNI',    'NACC',   'NACC'}, ...
    'transition', {'CU2MCI', 'MCI2Dem', 'CU2MCI', 'MCI2Dem'}, ...
    'core',       {1.08,     1.11,      2.296979, 1.899399}, ...
    'tail',       {6.6,      5.1,       11.632146, 6.937098}, ...
    'pTail',      {0.15,     0.15,      0.15,     0.15});

% Replicates and permutations
if strcmp(mode, 'real')
    C.N_PERM        = 1000;   % kinematics (bw_sensitivity)
    C.N_REPS_CAL    = 500;    % calibrated simulation (p_cal)
    C.N_PERM_CAL    = 500;
    C.N_REPS_CHECK  = 200;    % sensitivity checks of the calibration
    C.SIM_REPS_NULL = 500;    % generic simulation, type I error
    C.SIM_REPS_ALT  = 200;    % generic simulation, power
    C.SIM_PERM      = 199;
    C.SIM_BOOT      = 100;
    C.SHOW_CI       = true;   % subject bootstrap bands of the figures (300 resamples)
else                          % a few minutes in total
    C.N_PERM        = 199;
    C.N_REPS_CAL    = 50;
    C.N_PERM_CAL    = 99;
    C.N_REPS_CHECK  = 30;
    C.SIM_REPS_NULL = 20;
    C.SIM_REPS_ALT  = 10;
    C.SIM_PERM      = 99;
    C.SIM_BOOT      = 20;
    C.SHOW_CI       = false;
    C.ANALYSES      = {'A', 'V'};   % the A x V comparisons only in 'real'
end

onPath = strsplit(path, pathsep);
for d = {fullfile(C.root, 'src'), fullfile(C.root, 'experiments')}
    if ~any(strcmp(onPath, d{1})), addpath(d{1}); end
end
if ~exist(C.resultsDir, 'dir'), mkdir(C.resultsDir); end
end
