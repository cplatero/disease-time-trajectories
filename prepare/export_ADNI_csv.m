function export_ADNI_csv(srcDir, outDir)
%EXPORT_ADNI_CSV  Write the ADNI visit tables used by the paper as CSV.
%
%   export_ADNI_csv(srcDir, outDir)
%
%   Author-side script. It needs the original ADNI files, which cannot be
%   redistributed, and is not needed to run the experiments once the CSV
%   files exist. The CSV files must NOT be uploaded to GitHub
%   (data/*.csv is in .gitignore).
%
%   srcDir  folder from which the original scripts were run (the one with
%           ./aux_DPM/), so that '../data/' contains
%             ADNIMERGE_250321.mat, uwnpsychsum.xlsx,
%             PET_Amyloid_ADNI_260224.mat, ucd_wmh_data.xlsx
%   outDir  destination, normally <repo>/data
%
%   It calls prepare_ADNI_Amyloid (prepare/aux_DPM, unchanged) exactly as
%   plot_tau_AV did, for group 1 (CU -> MCI: CN, SMC and stable MCI) and
%   group 2 (MCI -> dementia: MCI and dementia), and writes
%       ADNI_CU2MCI.csv, ADNI_MCI2Dem.csv
%   with the columns of data/README.md:
%       ID            RID
%       AGE_visit     AGE + M/12 (years)
%       ClinGroup     1 stable, 2 converter, 3 next stage, from (Convert, stMCI)
%                     as in plot_tau_AV: CU2MCI (0,0) (1,0) (0,1);
%                     MCI2Dem (0,1) (1,1) (0,2); other rows NaN
%       MEM, EF       ADNI_MEM, ADNI_EF
%       CENTILOIDS    amyloid PET, Centiloids
%       WMH           TOTAL_WMH, ml
%       A_bl          AMYLOID_STATUS_bl > 0 (NaN without PET)
%       CENTILOIDS_bl, WMH_bl   first valid value of each subject
%       onset_age_clin  converters only: AGE + ConvertTime/12, the age at
%                     the first visit with the new diagnosis
%
%   Check after running: the subject counts printed here must match the
%   ADNI rows of the paper (Table S1 / Sec. 2).
%
%   Authorship: Carlos Platero (UPM), with the collaboration of Claude
%   (Anthropic).

if nargin < 2, error('export_ADNI_csv:args', 'Usage: export_ADNI_csv(srcDir, outDir)'); end
if ~exist(outDir, 'dir'), mkdir(outDir); end
if ~is_absolute(outDir), outDir = fullfile(pwd, outDir); end   % cd(srcDir) follows
here = pwd;  cleanupDir = onCleanup(@() cd(here));
addpath(fullfile(fileparts(mfilename('fullpath')), 'aux_DPM'));
cd(srcDir);

TR = {'CU2MCI', 'MCI2Dem'};
SV = {[0 0; 1 0; 0 1], [0 1; 1 1; 0 2]};          % (Convert, stMCI) of ClinGroup 1..3
for g = 1:2
    fprintf('\n[export ADNI] group %d (%s)\n', g, TR{g});
    A = prepare_ADNI_Amyloid(g, 'PET', true, false, false);

    ClinGroup = nan(height(A), 1);
    for k = 1:3
        ClinGroup(A.Convert == SV{g}(k, 1) & A.stMCI == SV{g}(k, 2)) = k;
    end
    onset = nan(height(A), 1);
    m = ClinGroup == 2;
    onset(m) = A.AGE(m) + A.ConvertTime(m) / 12;
    A_bl = double(A.AMYLOID_STATUS_bl > 0);
    A_bl(isnan(A.AMYLOID_STATUS_bl)) = NaN;

    T = table(double(A.RID), A.AGE + A.M / 12, ClinGroup, A.ADNI_MEM, A.ADNI_EF, ...
              A.CENTILOIDS, A.TOTAL_WMH, A_bl, A.CENTILOIDS_bl, A.WMH_bl, onset, ...
              'VariableNames', {'ID', 'AGE_visit', 'ClinGroup', 'MEM', 'EF', ...
              'CENTILOIDS', 'WMH', 'A_bl', 'CENTILOIDS_bl', 'WMH_bl', 'onset_age_clin'});
    f = fullfile(outDir, sprintf('ADNI_%s.csv', TR{g}));
    writetable(T, f);

    fprintf('  %s: %d visits, %d subjects\n', f, height(T), numel(unique(T.ID)));
    for k = 1:3
        fprintf('    ClinGroup %d: %d subjects\n', k, numel(unique(T.ID(T.ClinGroup == k))));
    end
    fprintf('    ClinGroup NaN: %d subjects | A_bl known: %d | WMH_bl known: %d\n', ...
        numel(unique(T.ID(isnan(T.ClinGroup)))), numel(unique(T.ID(~isnan(T.A_bl)))), ...
        numel(unique(T.ID(~isnan(T.WMH_bl)))));
end
end


function tf = is_absolute(p)
p = char(p);
tf = ~isempty(p) && (p(1) == '/' || p(1) == '\' || (numel(p) > 1 && p(2) == ':'));
end
