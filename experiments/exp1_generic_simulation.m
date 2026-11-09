function exp1_generic_simulation(mode)
%EXP1_GENERIC_SIMULATION  Simulation study with a known truth (Sec. 3.4).
%
%   exp1_generic_simulation('real')   % paper settings (several hours)
%   exp1_generic_simulation('demo')   % few replicates (minutes)
%
%   Paper items (check numbers against the final PDF):
%     Fig. 2        type I error by estimator          -> typeI/sim_typeI.png
%     Table 3       power of the permutation test       -> power/sim_power.tex
%     Tables S7-S8  type I error (by design and by bw)  -> typeI/sim_typeI.tex
%     bias and coverage (not in the paper)              -> power/sim_bias.tex,
%                                                          power/sim_coverage.tex
%   Table 2 of the paper lists the parameters of this simulation.
%   It uses no cohort data: it only needs src/sim_framework.m.
%
%   Design (as in the paper): type I error with n = 75 and 300 subjects
%   per subgroup, sigma_tau = 0 and 1.5 years, 15 % and 85 % of
%   single-visit subjects, bw = 2, 3 and 6; power for level and
%   divergence effects of 0.1-0.5 SD with sigma_tau = 1 and 15 %
%   single-visit subjects, bw = 3.
%
%   Authorship: Carlos Platero (UPM), with the collaboration of Claude
%   (Anthropic).

if nargin < 1, mode = 'real'; end
C = config(mode);
out = fullfile(C.resultsDir, 'exp1_generic_simulation');

sim_framework('design', 'typeI', ...
    'nGrid', [75 300], 'sigTauGrid', [0 1.5], 'pSingleGrid', [0.15 0.85], ...
    'bwGrid', [2 3 6], 'bw', C.BW, ...
    'nRepsNull', C.SIM_REPS_NULL, 'nPerm', C.SIM_PERM, ...
    'seed', C.SEED, 'outDir', fullfile(out, 'typeI'), 'bPlots', true);

sim_framework('design', 'power', ...
    'nGrid', [75 300], 'effGrid', [0.1 0.2 0.3 0.5], ...
    'sigTau', 1, 'pSingle', 0.15, 'bwGrid', C.BW, 'bw', C.BW, ...
    'nRepsAlt', C.SIM_REPS_ALT, 'nPerm', C.SIM_PERM, 'nBoot', C.SIM_BOOT, ...
    'seed', C.SEED, 'outDir', fullfile(out, 'power'), 'bPlots', true);
end
