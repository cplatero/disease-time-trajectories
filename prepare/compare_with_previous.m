function D = compare_with_previous(oldDir, newDir)
%COMPARE_WITH_PREVIOUS  Check that the CSV pipeline reproduces the old results.
%
%   D = compare_with_previous('C:\path\to\6b_DPM', 'results/real/exp3_bandwidth')
%
%   Compares bw_sensitivity_<C>_ll.csv and bw_sensitivity_<C>.csv of the
%   original scripts (oldDir) with those written by exp3_bandwidth
%   (newDir), comparison by comparison: subject counts, T_obs, p-value
%   and the start/end values of the curves. With the same data, seed and
%   number of permutations every difference must be 0 (up to rounding).
%   A difference in n0/n1 means that the CSV files do not select the same
%   subjects as the original loaders.
%
%   Only the comparisons present in both files are compared (the old
%   files also hold the A x V comparisons).
%
%   Authorship: Carlos Platero (UPM), with the collaboration of Claude
%   (Anthropic).

files = {'bw_sensitivity_ADNI_ll.csv', 'bw_sensitivity_NACC_ll.csv', ...
         'bw_sensitivity_ADNI.csv', 'bw_sensitivity_NACC.csv'};
cols  = {'n0', 'n1', 'T_obs', 'p', 'y_start0', 'y_end0', 'y_start1', 'y_end1'};
D = struct('file', {}, 'nCompared', {}, 'maxAbsDiff', {});
for f = 1:numel(files)
    fo = fullfile(oldDir, files{f});  fn = fullfile(newDir, files{f});
    if ~exist(fo, 'file') || ~exist(fn, 'file')
        fprintf('%-28s skipped (missing in one folder)\n', files{f});
        continue
    end
    O = readtable(fo, 'TextType', 'string');
    N = readtable(fn, 'TextType', 'string');
    keyO = make_key(O);  keyN = make_key(N);
    [~, io, in] = intersect(keyO, keyN);
    mx = zeros(1, numel(cols));
    for c = 1:numel(cols)
        d = abs(O.(cols{c})(io) - N.(cols{c})(in));
        d = d(~isnan(d));
        if ~isempty(d), mx(c) = max(d); end
    end
    D(end+1) = struct('file', files{f}, 'nCompared', numel(io), 'maxAbsDiff', mx); %#ok<AGROW>
    fprintf('%-28s %3d comparisons | max |diff|:', files{f}, numel(io));
    pairs = [cols; num2cell(mx)];
    fprintf(' %s=%.2g', pairs{:});
    fprintf('\n');
end
end


function k = make_key(T)
bw = arrayfun(@(b) sprintf('%g', b), T.bw, 'UniformOutput', false);
k = strcat(cellstr(T.transition), '|', cellstr(T.analysis), '|', cellstr(T.var), '|', bw);
end
