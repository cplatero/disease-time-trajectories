function plot_sim_cal(varargin)
%PLOT_SIM_CAL  Figure of the calibrated simulation (Fig. 3 of paper A).
%
%   plot_sim_cal()
%   plot_sim_cal('simFiles', {...}, 'pFiles', {...}, 'outFile', 'fig/fig_sim_cal.png')
%
%   Reads the results already computed and draws the figure. It does not
%   repeat any simulation or permutation test.
%     A. Type I error of the free permutation test in each real comparison,
%        for three estimators (local constant, local line, global
%        polynomial). The comparisons of both cohorts are ordered from
%        the smallest to the largest type I error of the local line.
%        Logarithmic scale. Dashed line: 0.05.
%     B. Permutation p-value against calibrated p-value, with local-linear
%        sw-LOESS. Dashed lines: 0.05. Solid line: identity.
%
%   Inputs (CSV files written by the analysis):
%     simFiles  sim_results.csv of sim_framework (design 'calibrated').
%               Columns used: scenario, estimator, reject, pCal.
%               scenario = cohort|transition|analysis|marker.
%     pFiles    bw_sensitivity_<cohort>_ll.csv of bw_sensitivity
%               (local line). Columns used: cohort, transition,
%               analysis, var, bw, p.
%
%   Options:
%     'bw'         bandwidth of the p-values in pFiles      [3]
%     'outFile'    PNG file ('' = do not save)              ['fig_sim_cal.png']
%     'dpi'        resolution of the PNG                    [300]
%     'showTitle'  overall title above the two panels       [false]
%
%   Authorship: Carlos Platero (UPM), with the collaboration of Claude
%   (Anthropic).

o = struct('simFiles', {{'sim_out_ADNI/sim_results.csv', 'sim_out_NACC/sim_results.csv'}}, ...
           'pFiles',   {{'bw_sensitivity_ADNI_ll.csv', 'bw_sensitivity_NACC_ll.csv'}}, ...
           'bw', 3, 'outFile', 'fig_sim_cal.png', 'dpi', 300, 'showTitle', false);
for k = 1:2:numel(varargin), o.(varargin{k}) = varargin{k+1}; end

EST  = {sprintf('swll_bw%g', o.bw), sprintf('sw_bw%g', o.bw), 'poly'};   % 1 = local line
ELAB = {'sw-LOESS, local line', 'sw-LOESS, local constant', 'Global polynomial'};
ECOL = [0.91 0.40 0.20; 0.18 0.47 0.85; 0.93 0.65 0.05];
EMRK = {'o', 'o', 's'};
DRAW = [2 1 3];            % drawing and legend order: constant, line, polynomial
COH  = {'ADNI', 'NACC'};   TRA = {'CU2MCI', 'MCI2Dem'};
TLAB = {'CU\rightarrowMCI', 'MCI\rightarrowDem'};
ANA  = {'A', 'V', 'AxV_A+', 'AxV_A-'};

% ---------- read the calibrated simulation (no duplicated scenarios)
S = struct('scenario', {{}}, 'estimator', {{}}, 'reject', [], 'pCal', []);
for f = 1:numel(o.simFiles)
    T = read_csv(o.simFiles{f});
    for r = 1:numel(T.scenario)
        dup = any(strcmp(S.scenario, T.scenario{r}) & strcmp(S.estimator, T.estimator{r}));
        if ~dup
            S.scenario{end+1} = T.scenario{r};  S.estimator{end+1} = T.estimator{r};
            S.reject(end+1)   = str2double(T.reject{r});
            S.pCal(end+1)     = str2double(T.pCal{r});
        end
    end
end

% ---------- read the permutation p-values (local line, bw = o.bw)
P = struct('key', {{}}, 'p', []);
for f = 1:numel(o.pFiles)
    T = read_csv(o.pFiles{f});
    bw = str2double(T.bw);
    for r = find(abs(bw - o.bw) < 1e-9)'
        P.key{end+1} = strjoin({T.cohort{r}, T.transition{r}, T.analysis{r}, T.var{r}}, '|');
        P.p(end+1)   = str2double(T.p{r});
    end
end

% ---------- order of the comparisons: cohort, transition, analysis, marker
keys = unique(S.scenario);
ord = zeros(numel(keys), 4);
for k = 1:numel(keys)
    q = strsplit(keys{k}, '|');
    ord(k, :) = [pos(COH, q{1}), pos(TRA, q{2}), pos(ANA, q{3}), k];
end
ord = sortrows(ord, [1 2 3 4]);
keys = keys(ord(:, 4));  grp = (ord(:, 1) - 1) * 2 + ord(:, 2);
nC = numel(keys);
fprintf('[plot_sim_cal] %d real comparisons.\n', nC);

typeI = nan(nC, numel(EST));  pcal = nan(nC, 1);  pperm = nan(nC, 1);
for k = 1:nC
    for e = 1:numel(EST)
        m = find(strcmp(S.scenario, keys{k}) & strcmp(S.estimator, EST{e}), 1);
        if ~isempty(m), typeI(k, e) = S.reject(m); end
        if e == 1 && ~isempty(m), pcal(k) = S.pCal(m); end
    end
    m = find(strcmp(P.key, keys{k}), 1);
    if ~isempty(m), pperm(k) = P.p(m); end
end
fprintf('[plot_sim_cal] Minimum calibrated p = %.3f. Comparisons without p: %d.\n', ...
    min(pcal), sum(isnan(pperm)));

% ---------- figure
fh = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 7.5]);
FS = 8;

% A. type I error per comparison, ordered by the local line (both cohorts)
ax1 = axes('Parent', fh, 'Position', [0.08 0.17 0.52 0.70]); hold(ax1, 'on');
yl = [0.008 1.2];
[~, io] = sortrows([typeI(:, 1), typeI(:, 2)]);       % ties broken by the local constant
h = gobjects(1, numel(EST) + 1);
for e = DRAW
    y = typeI(io, e);  y(y <= 0) = yl(1);           % zero rejections drawn at the axis
    h(e) = plot(ax1, 1:nC, y, EMRK{e}, 'MarkerSize', 4, 'LineStyle', 'none', ...
        'MarkerEdgeColor', 'w', 'MarkerFaceColor', ECOL(e, :));
end
h(end) = plot(ax1, [0.5 nC + 0.5], [0.05 0.05], 'k--', 'LineWidth', 0.8);
set(ax1, 'YScale', 'log', 'YLim', yl, 'XLim', [0 nC + 1], 'FontSize', FS, 'XTick', [], ...
    'YTick', [0.01 0.02 0.05 0.1 0.2 0.5 1], 'YTickLabel', {'0.01','0.02','0.05','0.1','0.2','0.5','1'}, ...
    'Box', 'off', 'TickDir', 'out');
grid(ax1, 'on'); set(ax1, 'XGrid', 'off', 'YMinorGrid', 'off');
ylabel(ax1, 'Type I error (\alpha = 0.05)');
xlabel(ax1, 'Real comparisons (ordered by the type I error of the local line)');
title(ax1, 'A. Type I error in each comparison', 'FontSize', FS + 1);
legend(h(DRAW), ELAB(DRAW), 'Location', 'northwest', 'FontSize', FS - 1, 'Box', 'off');

% B. permutation p against calibrated p (local line)
ax2 = axes('Parent', fh, 'Position', [0.70 0.17 0.28 0.70]); hold(ax2, 'on');
lim = [5e-4 1];
plot(ax2, lim, lim, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
plot(ax2, [0.05 0.05], lim, 'k--', 'LineWidth', 0.6, 'HandleVisibility', 'off');
plot(ax2, lim, [0.05 0.05], 'k--', 'LineWidth', 0.6, 'HandleVisibility', 'off');
CCOL = [0.12 0.47 0.82; 0.76 0.09 0.36];  TMRK = {'o', 's'};
hb = gobjects(1, 4);  lb = cell(1, 4);  n = 0;
for c = 1:2
    for t = 1:2
        m = grp == (c - 1) * 2 + t;
        n = n + 1;
        hb(n) = plot(ax2, max(pperm(m), lim(1)), max(pcal(m), lim(1)), TMRK{t}, 'MarkerSize', 4, ...
            'LineStyle', 'none', 'MarkerEdgeColor', CCOL(c, :), 'MarkerFaceColor', CCOL(c, :) * (t == 1) + [1 1 1] * (t == 2));
        lb{n} = sprintf('%s %s', COH{c}, TLAB{t});
    end
end
set(ax2, 'XScale', 'log', 'YScale', 'log', 'XLim', lim, 'YLim', lim, 'FontSize', FS, 'Box', 'on', ...
    'XTick', [0.001 0.01 0.05 0.2 1], 'XTickLabel', {'0.001','0.01','0.05','0.2','1'}, ...
    'YTick', [0.001 0.01 0.05 0.2 1], 'YTickLabel', {'0.001','0.01','0.05','0.2','1'});
axis(ax2, 'square'); grid(ax2, 'on'); set(ax2, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlabel(ax2, 'Permutation p-value');  ylabel(ax2, 'Calibrated p-value');
title(ax2, 'B. Local line: p vs. calibrated p', 'FontSize', FS + 1);
legend(hb, lb, 'Location', 'northwest', 'FontSize', FS - 1.5, 'Box', 'off');

if o.showTitle && exist('sgtitle', 'file')
    sgtitle(fh, sprintf('Calibrated simulation of the %d real comparisons (bw = %g years)', nC, o.bw), 'FontSize', FS + 2);
end
if ~isempty(o.outFile)
    set(fh, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 18 7.5]);
    print(fh, o.outFile, '-dpng', sprintf('-r%d', o.dpi));
    fprintf('[plot_sim_cal] Saved %s\n', o.outFile);
end
end


function i = pos(list, s)
% Position of s in list (numel(list) + 1 if absent).
i = find(strcmp(list, s), 1);
if isempty(i), i = numel(list) + 1; end
end


function T = read_csv(f)
% Minimal CSV reader (header + text cells). Works in MATLAB and Octave.
fid = fopen(f, 'r');
if fid < 0, error('plot_sim_cal:file', 'Cannot open %s', f); end
hdr = strsplit(strtrim(fgetl(fid)), ',');
C = textscan(fid, repmat('%s', 1, numel(hdr)), 'Delimiter', ',');
fclose(fid);
T = struct();
for k = 1:numel(hdr)
    T.(matlab.lang.makeValidName(hdr{k})) = C{k};
end
end
