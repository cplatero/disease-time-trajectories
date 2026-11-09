function cfg = load_cohort(cohort, varargin)
%LOAD_COHORT  Read one cohort (ADNI or NACC) from the CSV files of data/.
%
%   cfg = load_cohort('ADNI')
%   cfg = load_cohort('NACC', 'dataDir', 'data', 'tauDir', 'data/tau')
%
%   Single loader for every experiment. For each transition (CU2MCI and
%   MCI2Dem) it reads
%       <dataDir>/<COHORT>_<TRANSITION>.csv     one row per visit
%       <tauDir>/tau_<COHORT>_<TRANSITION>.csv  one row per subject
%   and places every visit on the disease-time axis,
%       t*_ij = age at visit - tau_hat_i,
%   where tau_hat_i is the mean of the 10 transition-age estimates of the
%   disease progression model (columns tau_1 ... tau_10).
%
%   Visit file columns (see data/README.md):
%       ID, AGE_visit, ClinGroup, MEM, EF, CENTILOIDS, WMH,
%       A_bl, CENTILOIDS_bl, WMH_bl, onset_age_clin
%   ClinGroup: 1 = stable, 2 = converter, 3 = already in the next stage
%   (NaN: not used in the clinical row, e.g. ADNI SMC). A_bl: baseline
%   amyloid status (1 = A+, 0 = A-, NaN = no PET). WMH_bl in ml.
%
%   Output: struct array cfg(1:2), one element per transition, with the
%   fields used by bw_sensitivity and plot_kinematics:
%       trans_code   'CU2MCI' | 'MCI2Dem'
%       event_label  'CU to MCI' | 'MCI to Demen.'
%       leg_row1     legend of the clinical row
%       select       all visits with a transition-age estimate
%       amyloid      visits with baseline amyloid status (A_Status_Basal)
%       wmh          visits with baseline WMH (A_Status_Basal = V+)
%   Every table has the columns SubjID and Time_to_Event.
%
%   This function replaces the cohort-specific loaders of the original
%   scripts (prepare_ADNI_Amyloid + join_estimation_tau for ADNI and
%   prepare_NACC_data for NACC). The CSV files are written by the scripts
%   of prepare/, so the selection of subjects and visits is the same.
%   As in the original loaders, ADNI keeps the visits of every subject
%   with an estimate, and NACC also drops visits without age.
%
%   Options
%     'dataDir'        folder of the visit files            ['data']
%     'tauDir'         folder of the transition-age files   ['data/tau']
%     'wmhThreshold'   V+ threshold on WMH_bl, ml           [10]
%     'verbose'        print the size of each table         [true]
%
%   Authorship: Carlos Platero (UPM), with the collaboration of Claude
%   (Anthropic).

p = inputParser;
addRequired(p, 'cohort', @(x) any(strcmpi(x, {'ADNI', 'NACC'})));
addParameter(p, 'dataDir', 'data');
addParameter(p, 'tauDir', fullfile('data', 'tau'));
addParameter(p, 'wmhThreshold', 10);
addParameter(p, 'verbose', true);
parse(p, cohort, varargin{:});
o = p.Results;
cohort = upper(char(o.cohort));

TR  = {'CU2MCI', 'MCI2Dem'};
EV  = {'CU to MCI', 'MCI to Demen.'};
if strcmp(cohort, 'ADNI'), third = 'Dementia'; else, third = 'sDem'; end
LEG = {{'sCU (stable)', 'pCU (converter)', 'sMCI', 'Global Trend'}, ...
       {'sMCI (stable)', 'pMCI (converter)', third, 'Global Trend'}};
REQ = {'ID', 'AGE_visit', 'ClinGroup', 'MEM', 'EF', 'CENTILOIDS', 'WMH', ...
       'A_bl', 'CENTILOIDS_bl', 'WMH_bl'};

cfg = struct('trans_code', {}, 'event_label', {}, 'leg_row1', {}, ...
             'select', {}, 'amyloid', {}, 'wmh', {});
for g = 1:2
    fv = fullfile(o.dataDir, sprintf('%s_%s.csv', cohort, TR{g}));
    ft = fullfile(o.tauDir,  sprintf('tau_%s_%s.csv', cohort, TR{g}));
    if ~exist(fv, 'file')
        error('load_cohort:file', ['Missing %s.\nWrite it with prepare/export_%s_csv.m, ' ...
              'or run the demo with run_all(''demo'').'], fv, cohort);
    end
    if ~exist(ft, 'file')
        error('load_cohort:file', 'Missing %s.\nWrite it with prepare/export_tau_csv.m.', ft);
    end

    D = readtable(fv);
    miss = setdiff(REQ, D.Properties.VariableNames);
    if ~isempty(miss)
        error('load_cohort:columns', '%s lacks the columns: %s', fv, strjoin(miss, ', '));
    end

    % Transition age: mean of the L = 10 model estimates of each subject.
    T = readtable(ft);
    vn = T.Properties.VariableNames;
    tauCols = strncmp(vn, 'tau_', 4);
    E = table2array(T(:, find(tauCols)));
    E(isinf(E)) = NaN;
    if size(E, 2) ~= 10
        warning('load_cohort:nEst', '%s has %d estimate columns (expected 10).', ft, size(E, 2));
    end
    tauMean = mean(E, 2, 'omitnan');

    [tf, loc] = ismember(D.ID, T.ID);
    tauVisit = nan(size(D.ID));
    tauVisit(tf) = tauMean(loc(tf));

    D.SubjID = D.ID;
    D.Time_to_Event = D.AGE_visit - tauVisit;
    D = D(~isnan(tauVisit), :);                    % subjects with an estimate
    if strcmp(cohort, 'NACC')
        D = D(~isnan(D.Time_to_Event), :);         % as prepare_NACC_data
    end

    A = D(~isnan(D.A_bl), :);
    A.A_Status_Basal = A.A_bl > 0;
    V = D(~isnan(D.WMH_bl), :);
    V.A_Status_Basal = V.WMH_bl > o.wmhThreshold;

    cfg(g).trans_code  = TR{g};
    cfg(g).event_label = EV{g};
    cfg(g).leg_row1    = LEG{g};
    cfg(g).select      = D;
    cfg(g).amyloid     = A;
    cfg(g).wmh         = V;

    if o.verbose
        fprintf('  %s %-7s select: %5d visits (%4d subjects) | A: %4d subj (A+ = %4d) | V: %4d subj (V+ = %4d)\n', ...
            cohort, TR{g}, height(D), numel(unique(D.SubjID)), ...
            numel(unique(A.SubjID)), numel(unique(A.SubjID(A.A_Status_Basal))), ...
            numel(unique(V.SubjID)), numel(unique(V.SubjID(V.A_Status_Basal))));
    end
end
end
