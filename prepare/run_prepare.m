% run_prepare.m
% Author-side step: converts the original ADNI and NACC tables and the
% DPM transition-age estimates into the CSV files read by the
% experiments. Run it once, from the root of the repository, after
% editing the two paths below. Nobody else needs to run it.
%
% The cohort CSV files (data/*.csv) are not versioned: the ADNI and NACC
% Data Use Agreements do not allow redistribution. The transition-age
% files (data/tau/*.csv) are not versioned either, until you confirm
% that both agreements allow publishing them (see data/README.md).
%
% Authorship: Carlos Platero (UPM), with the collaboration of Claude
% (Anthropic).

SRC_DIR  = 'C:\path\to\6b_DPM';                         % folder of the original scripts (../data reachable)
HIST_DIR = fullfile(SRC_DIR, '..', 'data', 'historia_natural');

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'prepare'));

export_tau_csv (HIST_DIR, fullfile(root, 'data', 'tau'));
export_NACC_csv(HIST_DIR, fullfile(root, 'data'));
export_ADNI_csv(SRC_DIR,  fullfile(root, 'data'));

% Check: the counts printed above must match the paper. Then run
%   run_all('real', 3)
%   compare_with_previous(SRC_DIR, fullfile(root, 'results', 'real', 'exp3_bandwidth'))
% Every difference must be 0: the CSV files then select the same subjects
% and visits as the original loaders. Finally run run_all('real') and
% compare results/real/ with the tables and figures of the paper.
