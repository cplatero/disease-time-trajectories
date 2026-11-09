function T = calib_tau_error(cohort, transition, varargin)
%CALIB_TAU_ERROR  Error model of the estimated transition age (Table 1).
%
%   T = calib_tau_error('NACC', 'CU2MCI')
%   T = calib_tau_error('ADNI', 'MCI2Dem', 'gap', 1, 'dataDir', 'data', 'tauDir', 'data/tau')
%
%   In the subjects who convert during follow-up (ClinGroup = 2), the
%   error is e_i = tau_hat_i - onset_age_clin_i, with tau_hat_i the mean
%   of the 10 DPM estimates. The clinical onset is the age at the first
%   visit with the new diagnosis, so the true onset lies in the preceding
%   visit interval g. Assuming it is uniform there, the clinical onset is
%   late by g/2 on average and adds a variance g^2/12. The error is
%   modelled as a mixture of two centred normals:
%     core (85 %): SD = robust SD of e (1.4826 MAD), corrected by g^2/12;
%     tail (15 %): SD such that 0.85 core^2 + 0.15 tail^2 equals the
%                  classical variance of e corrected by g^2/12.
%
%   Output T: struct with fields cohort, transition, n, bias, sd_classic,
%   sd_robust, core, tail, pTail, ready for the 'tauTable' of
%   sim_framework.
%
%   Inputs are the CSV files of data/ (load_cohort): the visit file gives
%   ClinGroup and onset_age_clin, the tau file the 10 estimates.
%
%   Options
%     'gap'      years between visits [1]
%     'pTail'    weight of the tail [0.15]
%     'dataDir'  ['data']    'tauDir'  ['data/tau']
%
%   Generalises the NACC version of 27 Oct 2026 (same formulas) to both
%   cohorts.
%
%   Authorship: Carlos Platero (UPM), with the collaboration of Claude
%   (Anthropic).

p = inputParser;
addParameter(p, 'gap', 1);
addParameter(p, 'pTail', 0.15);
addParameter(p, 'dataDir', 'data');
addParameter(p, 'tauDir', fullfile('data', 'tau'));
parse(p, varargin{:});
o = p.Results;
cohort = upper(char(cohort));  transition = char(transition);

D = readtable(fullfile(o.dataDir, sprintf('%s_%s.csv', cohort, transition)));
[ids, first] = unique(D.ID, 'stable');            % one row per subject
grp   = D.ClinGroup(first);
onset = D.onset_age_clin(first);

A  = readtable(fullfile(o.tauDir, sprintf('tau_%s_%s.csv', cohort, transition)));
vn = A.Properties.VariableNames;
E  = table2array(A(:, find(strncmp(vn, 'tau_', 4))));
E(isinf(E)) = NaN;
tau_mean = mean(E, 2, 'omitnan');

[~, it, id] = intersect(ids, A.ID, 'stable');
m = grp(it) == 2;
e = tau_mean(id(m)) - onset(it(m));
e = e(~isnan(e));
n = numel(e);
if n < 10
    error('calib_tau_error:n', '%s %s: only %d converters with tau_hat and clinical onset.', cohort, transition, n);
end

g2   = o.gap^2 / 12;
sd_c = std(e);
sd_r = 1.4826 * median(abs(e - median(e)));
core = sqrt(max(sd_r^2 - g2, 0));
var_c = max(sd_c^2 - g2, 0);
tail = sqrt(max((var_c - (1 - o.pTail) * core^2) / o.pTail, core^2));

T = struct('cohort', cohort, 'transition', transition, 'n', n, ...
           'bias', mean(e) + o.gap / 2, 'sd_classic', sd_c, 'sd_robust', sd_r, ...
           'core', core, 'tail', tail, 'pTail', o.pTail);

fprintf('\n%s %s: %d converters with tau_hat and clinical onset\n', cohort, transition, n);
fprintf('  bias (tau_hat - onset + g/2) = %+.2f years\n', T.bias);
fprintf('  SD classic = %.2f | SD robust (1.4826 MAD) = %.2f | IQR/1.349 = %.2f\n', ...
        sd_c, sd_r, iqr(e) / 1.349);
fprintf('  mixture: core SD = %.2f (%.0f %%), tail SD = %.2f (%.0f %%)\n', ...
        core, 100 * (1 - o.pTail), tail, 100 * o.pTail);
end
