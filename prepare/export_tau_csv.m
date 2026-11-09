function export_tau_csv(histDir, outDir)
%EXPORT_TAU_CSV  Write the transition-age estimates of the DPMs as CSV.
%
%   export_tau_csv(histDir, outDir)
%
%   histDir  folder with the age_onset_*.mat files
%            (originally ../data/historia_natural/)
%   outDir   destination, normally <repo>/data/tau
%
%   Each age_onset matrix has the subject ID in column 1 and the L = 10
%   Leaspy estimates of the transition age in columns 2:11. Inf is
%   written as NaN, and subjects without any estimate are dropped (as in
%   load_tau_estimations / prepare_NACC_data of the original scripts).
%
%       age_onset_leaspy578_preAD.mat         -> tau_ADNI_CU2MCI.csv
%       age_onset_leaspy_MCI2AD_7_8.mat       -> tau_ADNI_MCI2Dem.csv
%       age_onset_NACC_CU2MCI_2_3_4_13.mat    -> tau_NACC_CU2MCI.csv
%       age_onset_NACC_MCI2Dem_1_2_3_13.mat   -> tau_NACC_MCI2Dem.csv
%
%   Before publishing these files, check that the ADNI and NACC Data Use
%   Agreements allow it: they contain one derived value per subject ID.
%
%   Authorship: Carlos Platero (UPM), with the collaboration of Claude
%   (Anthropic).

if nargin < 2, error('export_tau_csv:args', 'Usage: export_tau_csv(histDir, outDir)'); end
if ~exist(outDir, 'dir'), mkdir(outDir); end

IN  = {'age_onset_leaspy578_preAD.mat', 'age_onset_leaspy_MCI2AD_7_8.mat', ...
       'age_onset_NACC_CU2MCI_2_3_4_13.mat', 'age_onset_NACC_MCI2Dem_1_2_3_13.mat'};
OUT = {'tau_ADNI_CU2MCI.csv', 'tau_ADNI_MCI2Dem.csv', ...
       'tau_NACC_CU2MCI.csv', 'tau_NACC_MCI2Dem.csv'};
for k = 1:numel(IN)
    S = load(fullfile(histDir, IN{k}), 'age_onset');
    A = double(S.age_onset);
    A(isinf(A)) = NaN;
    if size(A, 2) ~= 11
        error('export_tau_csv:layout', '%s has %d columns (expected ID + 10 estimates).', IN{k}, size(A, 2));
    end
    A = A(~all(isnan(A(:, 2:end)), 2), :);
    names = [{'ID'}, arrayfun(@(j) sprintf('tau_%d', j), 1:10, 'UniformOutput', false)];
    T = array2table(A, 'VariableNames', names);
    writetable(T, fullfile(outDir, OUT{k}));
    fprintf('[export tau] %-38s -> %-22s %5d subjects\n', IN{k}, OUT{k}, height(T));
end
end
