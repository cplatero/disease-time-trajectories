function export_NACC_csv(histDir, outDir, clThreshold)
%EXPORT_NACC_CSV  Write the NACC visit tables used by the paper as CSV.
%
%   export_NACC_csv(histDir, outDir)
%   export_NACC_csv(histDir, outDir, 20)
%
%   Author-side script. NACC data cannot be redistributed: the CSV files
%   must NOT be uploaded to GitHub (data/*.csv is in .gitignore).
%
%   histDir      folder with nacc_long_NACC_<TRANSITION>_<combo>.mat
%                (originally ../data/historia_natural/)
%   outDir       destination, normally <repo>/data
%   clThreshold  A+ threshold on CENTILOIDS_bl, Centiloids [20]
%
%   It reads the same files as prepare_NACC_data of plot_tau_AV
%       nacc_long_NACC_CU2MCI_2_3_4_13.mat
%       nacc_long_NACC_MCI2Dem_1_2_3_13.mat
%   and writes NACC_CU2MCI.csv and NACC_MCI2Dem.csv with the columns of
%   data/README.md. ClinGroup comes from (group, CONVERT) as in
%   plot_tau_AV: (0,0) stable, (0,1) converter, (1,0) next stage.
%   A_bl = CENTILOIDS_bl > clThreshold (NaN without Centiloids), the
%   rule of prepare_NACC_data.
%
%   Authorship: Carlos Platero (UPM), with the collaboration of Claude
%   (Anthropic).

if nargin < 2, error('export_NACC_csv:args', 'Usage: export_NACC_csv(histDir, outDir)'); end
if nargin < 3 || isempty(clThreshold), clThreshold = 20; end
if ~exist(outDir, 'dir'), mkdir(outDir); end

TR    = {'CU2MCI', 'MCI2Dem'};
COMBO = {'2_3_4_13', '1_2_3_13'};
SV    = [0 0; 0 1; 1 0];                      % (group, CONVERT) of ClinGroup 1..3
for g = 1:2
    f_in = fullfile(histDir, sprintf('nacc_long_NACC_%s_%s.mat', TR{g}, COMBO{g}));
    S = load(f_in, 'nacc_long');
    nl = S.nacc_long;
    col = @(name) double(nl.(name)(:));

    grp = col('group');  conv = col('CONVERT');
    ClinGroup = nan(size(grp));
    for k = 1:3
        ClinGroup(grp == SV(k, 1) & conv == SV(k, 2)) = k;
    end
    CLbl = col('CENTILOIDS_bl');
    A_bl = double(CLbl > clThreshold);
    A_bl(isnan(CLbl)) = NaN;
    onset = col('onset_age_clin');
    onset(ClinGroup ~= 2) = NaN;              % only converters are used

    T = table(col('ID'), col('AGE_visit'), ClinGroup, col('MEM'), col('EF'), ...
              col('CENTILOIDS'), col('WMH'), A_bl, CLbl, col('WMH_bl'), onset, ...
              'VariableNames', {'ID', 'AGE_visit', 'ClinGroup', 'MEM', 'EF', ...
              'CENTILOIDS', 'WMH', 'A_bl', 'CENTILOIDS_bl', 'WMH_bl', 'onset_age_clin'});
    f = fullfile(outDir, sprintf('NACC_%s.csv', TR{g}));
    writetable(T, f);

    fprintf('[export NACC] %s: %d visits, %d subjects\n', f, height(T), numel(unique(T.ID)));
    for k = 1:3
        fprintf('    ClinGroup %d: %d subjects\n', k, numel(unique(T.ID(T.ClinGroup == k))));
    end
    fprintf('    A_bl known: %d | WMH_bl known: %d\n', ...
        numel(unique(T.ID(~isnan(T.A_bl)))), numel(unique(T.ID(~isnan(T.WMH_bl)))));
end
end
