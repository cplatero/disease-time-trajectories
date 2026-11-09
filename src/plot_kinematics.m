% plot_kinematics.m
%
% Kinematics of the markers on the disease-time axis (Figs. 4-5 and
% Tables S11-S12 of the paper). One figure per cohort and transition:
%   row 1  clinical groups (stable, converter, next stage) and the global
%          curve of the transition;
%   row 2  baseline amyloid status, A- vs A+;
%   row 3  baseline vascular burden, V- vs V+;
%   columns MEM, EF, Centiloids and WMH.
% Every curve is a subject-weighted local-linear sw-LOESS curve
% (bw = 3 years) drawn on its coverage range, with a 95 % subject
% bootstrap band. Each pair of subgroup curves is compared with the
% subject-level permutation test on their common range.
%
% USAGE
%   plot_kinematics('ADNI', 2, true, 'pFile', 'bw_sensitivity_ADNI_ll.csv', ...
%       'pcalFile', 'sim_out_ADNI/sim_results.csv', 'saveDir', 'fig', ...
%       'dataDir', 'data', 'tauDir', 'data/tau', 'outDir', 'results')
%
% INPUTS
%   cohort       'ADNI' | 'NACC'
%   poly_degree  degree of the polynomial curves; only used with
%                'allLoess', false (curves of the original scripts) [2]
%   show_CI      draw the bootstrap bands [true]
% NAME-VALUE OPTIONS
%   'estimator'  'linear' (local-linear sw-LOESS, paper) | 'const' (local constant)
%   'pFile'      bw_sensitivity_<cohort>_ll.csv: permutation p-values of the
%                paper (same permutations as the bandwidth analysis).
%                Comparisons not found there (the marker that defines each
%                stratification) are tested here, 1000 permutations.
%   'pcalFile'   sim_results.csv of the calibrated simulation (p_cal)
%   'saveDir'    folder of the figures, HN_<cohort>_<CU2MCI|MCI2AD>_<ver>.png
%   'figVersion' suffix of the figure names ['v3']
%   'allLoess'   true: sw-LOESS in every subplot (paper) [true]
%   'dataDir', 'tauDir', 'wmhThreshold'   passed to load_cohort
%   'outDir'     folder of the LaTeX table metrics_table_<cohort>_ll_v3.tex ['.']
%
% Derived from plot_tau_AV.m (v3, 29 Oct 2026). The only changes are the
% data source (load_cohort, CSV files) and the language of the table.
%
% Authorship: Carlos Platero (UPM), with the collaboration of Claude
% (Anthropic).

function plot_kinematics(cohort, poly_degree, show_CI, varargin)

    if nargin < 1 || isempty(cohort)
        error('plot_kinematics:cohort', 'cohort is required: ''ADNI'' or ''NACC''.');
    end
    if nargin < 2 || isempty(poly_degree), poly_degree = 2;    end
    if nargin < 3 || isempty(show_CI),     show_CI     = true;  end
    opt = parse_options(varargin{:});

    cohort = upper(string(cohort));
    if ~any(strcmp(char(cohort), {'ADNI', 'NACC'}))
        error('plot_kinematics:cohort', 'cohort must be ''ADNI'' or ''NACC''.');
    end
    % CENTILOIDS and WMH (columns 3-4) are drawn as a point cloud
    % (plot_scatter_local_h), not as spaghetti lines, regardless of
    % cohort: most subjects have a single valid CL/WMH measurement, so
    % plot_lines_local_h (which needs >=2 points per subject to draw
    % anything) silently drops them, and the plot ends up showing only
    % the population trend, no individual data. This was originally
    % NACC-only (single-timepoint imaging is the norm there); it also
    % happened in ADNI once the doubly-stratified A x V subgroups made
    % the per-subject visit count for CL/WMH sparse there too.

    %% ================================================================
    %% GLOBAL PARAMETERS (identical to both original scripts)
    %% ================================================================
    LOESS_BW       = 3.0;
    DENSITY_TRIM_P = 10;
    N_BOOT         = 300;
    N_QUERY        = 150;
    N_PERM         = 1000;
    T_RANGE        = [-10, 5];
    VARS           = {'MEM', 'EF', 'CENTILOIDS', 'WMH'};
    SUBGRP_RANGE_P = 10;
    rng(2026,'twister');
    ROW_LABEL_FONTSIZE = 11;   % font size of the row label (ylabel,
                                % 1st column: 'Clinical Progression',
                                % 'Basal Amyloid Status', ...)

    USE_LOESS = false(3, 4);
    USE_LOESS(2, 1:2) = true;
    USE_LOESS(2, 4)   = true;
    USE_LOESS(3, 1:2) = true;
    USE_LOESS(3, 3)   = true;
    if opt.allLoess
        USE_LOESS = true(3, 4);    % v3: sw-LOESS in every subplot
    end

    % Subgroup smoother: local linear (default) or local constant.
    if strcmp(opt.estimator, 'linear')
        SMOOTHER = @sw_loess_ll_eval;   est_tag = 'local-linear sw-LOESS';
    else
        SMOOTHER = @sw_loess_eval;      est_tag = 'local-constant sw-LOESS';
    end
    PT  = read_p_table(opt.pFile, opt.estimator, LOESS_BW);
    PCT = read_pcal_table(opt.pcalFile, opt.estimator, LOESS_BW);
    fprintf('[AV] Estimator: %s | p-values: %s | p_cal: %s\n', est_tag, ...
        ternary(isempty(PT), 'computed here', opt.pFile), ...
        ternary(isempty(PCT), 'none', opt.pcalFile));

    %% ================================================================
    %% STEP 1 -- Load data + global Y limits
    %% ================================================================
    fprintf('\n[AV] Loading %s data...\n', cohort);

    cfg = load_cohort(char(cohort), 'dataDir', opt.dataDir, 'tauDir', opt.tauDir, ...
                      'wmhThreshold', opt.wmhThreshold);

    all_vals_for_ylim = cell(4, 1);
    for g = 1:numel(cfg)
        for i = 1:4
            vn = VARS{i};
            v = [cfg(g).select.(vn); cfg(g).amyloid.(vn); cfg(g).wmh.(vn)];
            all_vals_for_ylim{i} = [all_vals_for_ylim{i}; v(~isnan(v))];
        end
    end
    y_lims = zeros(4, 2);
    for i = 1:4
        v = all_vals_for_ylim{i};
        if isempty(v), y_lims(i,:) = [0 1]; continue; end
        lo = prctile(v, 2);  hi = prctile(v, 98);
        mg = 0.05*(hi-lo);
        y_lims(i,:) = [lo-mg, hi+mg];
    end

    %% ================================================================
    %% STEP 2 -- Figures + metrics collection
    %% ================================================================
    metrics = struct('group',{}, 'row',{}, 'col',{}, 'var',{}, ...
                     'subgrp',{}, 't_start',{}, 't_end',{}, ...
                     'y_start',{}, 'y_end',{}, 'slope',{}, ...
                     'pval',{}, 'n_subj',{}, 'pcal',{});

    for g = 1:numel(cfg)
        ev         = cfg(g).event_label;
        select_tbl = cfg(g).select;
        amyloid_tbl= cfg(g).amyloid;
        wmh_tbl    = cfg(g).wmh;
        leg_row1   = cfg(g).leg_row1;
        tcode      = cfg(g).trans_code;

        fprintf('\n[AV] Drawing figure %d (%s, %s)...\n', g, cohort, ev);
        fprintf('  select: %d visits (%d subjects) | amyloid: %d subj (A+=%d) | wmh: %d subj (V+=%d)\n', ...
            height(select_tbl), numel(unique(select_tbl.SubjID)), ...
            numel(unique(amyloid_tbl.SubjID)), sum(amyloid_tbl.A_Status_Basal), ...
            numel(unique(wmh_tbl.SubjID)), sum(wmh_tbl.A_Status_Basal));

        row_labels = { ['Clinical Progression  (' ev ')'], ...
                       ['Basal Amyloid Status  (' ev ')'], ...
                       ['Basal Vascular Burden (' ev ')'] };

        figure('Color','w','Name',[char(cohort) ' Natural History: ' ev], ...
               'Units','normalized','Position',[0.02 0.05 0.96 0.88]);
        ax = gobjects(3,4);

        %% -- ROW 1: clinical groups -----------------------------------
        for i = 1:4
            vn = VARS{i};
            ax(1,i) = subplot(3,4,i);
            hold on; grid on; box on;
            xlim(T_RANGE);
            if y_lims(i,1) < y_lims(i,2), ylim(y_lims(i,:)); end

            plot_fn = @plot_lines_local_h;
            if i > 2, plot_fn = @plot_scatter_local_h; end

            h1 = plot_fn(select_tbl(select_tbl.ClinGroup==1,:), vn, [0.00 0.80 0.00]);
            h2 = plot_fn(select_tbl(select_tbl.ClinGroup==2,:), vn, [0.00 0.40 1.00]);
            h3 = plot_fn(select_tbl(select_tbl.ClinGroup==3,:), vn, [0.80 0.00 0.00]);

            if USE_LOESS(1,i)
                % Global curve: same estimator and coverage rule as the
                % subgroup curves, applied to the whole transition.
                trG = subgrp_range(select_tbl, SUBGRP_RANGE_P, T_RANGE);
                [h4, y_plot, x_plot] = plot_trend_loess_sw_out(select_tbl, vn, ...
                        [0 0 0], 4, LOESS_BW, trG, show_CI, N_BOOT, N_QUERY, SMOOTHER);
            else
                [h4, y_plot, x_plot] = plot_trend_poly_out(select_tbl, vn, ...
                        [0 0 0], 4, poly_degree, DENSITY_TRIM_P, T_RANGE, ...
                        show_CI, N_BOOT, N_QUERY);
            end

            xline(0,'--k','LineWidth',1.5,'HandleVisibility','off');
            title(vn,'FontWeight','bold');
            xlabel('Years to Event');
            if i==1, ylabel(row_labels{1},'FontWeight','bold','FontSize',ROW_LABEL_FONTSIZE); end
            if i==4
                % This subplot (i==4, WMH) uses plot_scatter_local_h for
                % h1-h3, so the legend describes markers, not lines.
                leg4 = leg_row1;
                leg4{1} = [leg4{1} ' (pts)'];
                leg4{2} = [leg4{2} ' (pts)'];
                leg4{3} = [leg4{3} ' (pts)'];
                build_legend([h1 h2 h3 h4], leg4);
            end

            if ~isempty(y_plot)
                m = make_metric(g, 1, i, vn, 'Global', x_plot, y_plot, NaN, ...
                        numel(unique(select_tbl.SubjID( ...
                            ~isnan(select_tbl.(vn)) & ~isnan(select_tbl.Time_to_Event)))), NaN);
                metrics(end+1) = m; %#ok<AGROW>
            end
        end

        %% -- ROW 2: baseline amyloid status -----------------------------
        grp0_A = amyloid_tbl(amyloid_tbl.A_Status_Basal==0,:);
        grp1_A = amyloid_tbl(amyloid_tbl.A_Status_Basal==1,:);

        for i = 1:4
            vn = VARS{i};
            ax(2,i) = subplot(3,4,4+i);
            hold on; grid on; box on;
            xlim(T_RANGE);
            if y_lims(i,1) < y_lims(i,2), ylim(y_lims(i,:)); end

            plot_fn = @plot_lines_local_h;
            if i > 2, plot_fn = @plot_scatter_local_h; end
            h1 = plot_fn(grp0_A, vn, [0.00 0.80 0.00]);
            h2 = plot_fn(grp1_A, vn, [1.00 0.00 0.00]);

            tr0 = subgrp_range(grp0_A, SUBGRP_RANGE_P, T_RANGE);
            tr1 = subgrp_range(grp1_A, SUBGRP_RANGE_P, T_RANGE);

            if USE_LOESS(2,i)
                [h3, yp0, xp0] = plot_trend_loess_sw_out(grp0_A, vn, ...
                        [0.00 0.70 0.00], 4, LOESS_BW, tr0, show_CI, N_BOOT, N_QUERY, SMOOTHER);
                [h4, yp1, xp1] = plot_trend_loess_sw_out(grp1_A, vn, ...
                        [1.00 0.00 0.00], 4, LOESS_BW, tr1, show_CI, N_BOOT, N_QUERY, SMOOTHER);
            else
                [h3, yp0, xp0] = plot_trend_poly_out(grp0_A, vn, ...
                        [0.00 0.70 0.00], 4, poly_degree, DENSITY_TRIM_P, tr0, ...
                        show_CI, N_BOOT, N_QUERY);
                [h4, yp1, xp1] = plot_trend_poly_out(grp1_A, vn, ...
                        [1.00 0.00 0.00], 4, poly_degree, DENSITY_TRIM_P, tr1, ...
                        show_CI, N_BOOT, N_QUERY);
            end

            xline(0,'--k','LineWidth',1.5,'HandleVisibility','off');
            xlabel('Years to Event');
            if i==1, ylabel(row_labels{2},'FontWeight','bold','FontSize',ROW_LABEL_FONTSIZE); end
            if i==4
                build_legend([h1 h2 h3 h4], ...
                    {'A- Indiv.','A+ Indiv.','A- Trend','A+ Trend'});
            end

            pval = lookup_p(PT, tcode, 'A', vn);
            if isnan(pval)
                pval = curve_perm_test(grp0_A, grp1_A, vn, USE_LOESS(2,i), ...
                           LOESS_BW, tr0, tr1, poly_degree, DENSITY_TRIM_P, N_QUERY, N_PERM, SMOOTHER);
            end
            pcal = lookup_pcal(PCT, cohort, tcode, 'A', vn);
            for si = 0:1
                if si==0, yp=yp0; xp=xp0; grp=grp0_A; lbl='A-';
                else,      yp=yp1; xp=xp1; grp=grp1_A; lbl='A+'; end
                if isempty(yp), continue; end
                m = make_metric(g, 2, i, vn, lbl, xp, yp, pval, ...
                        numel(unique(grp.SubjID( ...
                            ~isnan(grp.(vn)) & ~isnan(grp.Time_to_Event)))), pcal);
                metrics(end+1) = m; %#ok<AGROW>
            end
        end

        %% -- ROW 3: baseline vascular burden -----------------------------
        grp0_V = wmh_tbl(wmh_tbl.A_Status_Basal==0,:);   % V-
        grp1_V = wmh_tbl(wmh_tbl.A_Status_Basal==1,:);   % V+

        for i = 1:4
            vn = VARS{i};
            ax(3,i) = subplot(3,4,8+i);
            hold on; grid on; box on;
            xlim(T_RANGE);
            if y_lims(i,1) < y_lims(i,2), ylim(y_lims(i,:)); end

            plot_fn = @plot_lines_local_h;
            if i > 2, plot_fn = @plot_scatter_local_h; end
            h1 = plot_fn(grp0_V, vn, [0.00 0.80 0.00]);
            h2 = plot_fn(grp1_V, vn, [1.00 0.00 0.00]);

            tr0 = subgrp_range(grp0_V, SUBGRP_RANGE_P, T_RANGE);
            tr1 = subgrp_range(grp1_V, SUBGRP_RANGE_P, T_RANGE);

            if USE_LOESS(3,i)
                [h3, yp0, xp0] = plot_trend_loess_sw_out(grp0_V, vn, ...
                        [0.00 0.70 0.00], 4, LOESS_BW, tr0, show_CI, N_BOOT, N_QUERY, SMOOTHER);
                [h4, yp1, xp1] = plot_trend_loess_sw_out(grp1_V, vn, ...
                        [1.00 0.00 0.00], 4, LOESS_BW, tr1, show_CI, N_BOOT, N_QUERY, SMOOTHER);
            else
                [h3, yp0, xp0] = plot_trend_poly_out(grp0_V, vn, ...
                        [0.00 0.70 0.00], 4, poly_degree, DENSITY_TRIM_P, tr0, ...
                        show_CI, N_BOOT, N_QUERY);
                [h4, yp1, xp1] = plot_trend_poly_out(grp1_V, vn, ...
                        [1.00 0.00 0.00], 4, poly_degree, DENSITY_TRIM_P, tr1, ...
                        show_CI, N_BOOT, N_QUERY);
            end

            xline(0,'--k','LineWidth',1.5,'HandleVisibility','off');
            xlabel('Years to Event');
            if i==1, ylabel(row_labels{3},'FontWeight','bold','FontSize',ROW_LABEL_FONTSIZE); end
            if i==4
                build_legend([h1 h2 h3 h4], ...
                    {'V- Indiv.','V+ Indiv.','V- Trend','V+ Trend'});
            end

            pval = lookup_p(PT, tcode, 'V', vn);
            if isnan(pval)
                pval = curve_perm_test(grp0_V, grp1_V, vn, USE_LOESS(3,i), ...
                           LOESS_BW, tr0, tr1, poly_degree, DENSITY_TRIM_P, N_QUERY, N_PERM, SMOOTHER);
            end
            pcal = lookup_pcal(PCT, cohort, tcode, 'V', vn);
            for si = 0:1
                if si==0, yp=yp0; xp=xp0; grp=grp0_V; lbl='V-';
                else,      yp=yp1; xp=xp1; grp=grp1_V; lbl='V+'; end
                if isempty(yp), continue; end
                m = make_metric(g, 3, i, vn, lbl, xp, yp, pval, ...
                        numel(unique(grp.SubjID( ...
                            ~isnan(grp.(vn)) & ~isnan(grp.Time_to_Event)))), pcal);
                metrics(end+1) = m; %#ok<AGROW>
            end
        end

        for i=1:4, linkaxes(ax(:,i),'y'); end

        if opt.allLoess
            ttl = sprintf('%s Natural History: %s  [%s, all curves | bw=%.0f | range_p=%d | CI=%s]', ...
                char(cohort), ev, est_tag, LOESS_BW, SUBGRP_RANGE_P, mat2str(show_CI));
        else
            ttl = sprintf('%s Natural History: %s  [%s | poly=%d | bw=%.0f | range_p=%d | CI=%s]', ...
                char(cohort), ev, est_tag, poly_degree, LOESS_BW, SUBGRP_RANGE_P, mat2str(show_CI));
        end
        sgtitle(ttl, 'FontSize',13,'FontWeight','bold');

        if ~isempty(opt.saveDir)
            fig_codes = {'CU2MCI', 'MCI2AD'};          % names used in paper A
            save_figure(gcf, fullfile(opt.saveDir, sprintf('HN_%s_%s_%s.png', ...
                char(cohort), fig_codes{g}, opt.figVersion)));
        end
    end

    %% ================================================================
    %% STEP 3 -- Build and save LaTeX table
    %% ================================================================
    if ~exist(opt.outDir, 'dir'), mkdir(opt.outDir); end
    fname = fullfile(opt.outDir, sprintf('metrics_table_%s%s%s.tex', char(cohort), ...
        ternary(strcmp(opt.estimator,'linear'), '_ll', ''), ternary(opt.allLoess, '_v3', '')));
    write_latex_table(metrics, fname, cohort, LOESS_BW, est_tag, ~isempty(PCT), opt.allLoess);
    fprintf('\n[AV] LaTeX table saved to %s\n', fname);
    fprintf('[AV] All done.\n');
end


%% ================================================================
%%  METRICS: permutation test on mean |curve_A - curve_B|
%%  (Eqs. 7-8 of the manuscript; +1 correction added, see header;
%%  uses the generic
%%  'SubjID' column set by load_cohort)
%% ================================================================
function pval = curve_perm_test(grp0, grp1, varName, use_loess, ...
                                 bw, tr0, tr1, poly_deg, trim_p, n_q, n_perm, smoother)
    t_lo = max(tr0(1), tr1(1));
    t_hi = min(tr0(2), tr1(2));
    if t_hi <= t_lo
        pval = NaN;
        return
    end
    x_common = linspace(t_lo, t_hi, n_q);

    y0_obs = eval_curve(grp0, varName, use_loess, bw, tr0, poly_deg, trim_p, x_common, smoother);
    y1_obs = eval_curve(grp1, varName, use_loess, bw, tr1, poly_deg, trim_p, x_common, smoother);
    ok = ~isnan(y0_obs) & ~isnan(y1_obs);
    if sum(ok) < 5, pval = NaN; return; end
    T_obs = mean(abs(y1_obs(ok) - y0_obs(ok)));

    all_data = [grp0; grp1];
    n0 = numel(unique(grp0.SubjID));
    n1 = numel(unique(grp1.SubjID));
    u_all = unique(all_data.SubjID);
    n_tot = numel(u_all);

    T_perm = zeros(n_perm, 1);
    for p = 1:n_perm
        idx    = randperm(n_tot);
        ids_p0 = u_all(idx(1:n0));
        ids_p1 = u_all(idx(n0+1:n0+n1));
        d0 = all_data(ismember(all_data.SubjID, ids_p0), :);
        d1 = all_data(ismember(all_data.SubjID, ids_p1), :);
        t_range_all = [min(all_data.Time_to_Event) max(all_data.Time_to_Event)];
        tr0p = subgrp_range(d0, 10, t_range_all);
        tr1p = subgrp_range(d1, 10, t_range_all);
        y0p = eval_curve(d0, varName, use_loess, bw, tr0p, poly_deg, trim_p, x_common, smoother);
        y1p = eval_curve(d1, varName, use_loess, bw, tr1p, poly_deg, trim_p, x_common, smoother);
        ok2 = ~isnan(y0p) & ~isnan(y1p);
        if sum(ok2) < 3
            T_perm(p) = 0;
        else
            T_perm(p) = mean(abs(y1p(ok2) - y0p(ok2)));
        end
    end
    % Eq. (8) of the manuscript (Phipson & Smyth, 2010): +1 in numerator
    % and denominator, so the p-value is never 0 (minimum 1/(1+n_perm)).
    pval = (1 + sum(T_perm >= T_obs)) / (1 + n_perm);
end


function y = eval_curve(data, varName, use_loess, bw, tr, poly_deg, trim_p, x_q, smoother)
    validMask = ~isnan(data.(varName)) & ~isnan(data.Time_to_Event);
    t_all  = data.Time_to_Event(validMask);
    y_all  = data.(varName)(validMask);
    id_all = data.SubjID(validMask);
    u_ids  = unique(id_all);
    if numel(u_ids) < 3
        y = nan(size(x_q));
        return
    end
    if use_loess
        y = smoother(t_all, y_all, id_all, u_ids, x_q, bw);
    else
        t_lo = max(tr(1), prctile(t_all, trim_p));
        t_hi = min(tr(2), prctile(t_all, 100-trim_p));
        if t_hi <= t_lo, y = nan(size(x_q)); return; end
        mask = t_all >= t_lo & t_all <= t_hi;
        if sum(mask) < poly_deg+2, y = nan(size(x_q)); return; end
        p = polyfit(t_all(mask), y_all(mask), poly_deg);
        y = polyval(p, x_q);
        y(x_q < t_lo | x_q > t_hi) = NaN;
    end
end


%% ================================================================
%%  LaTeX TABLE WRITER (format of tab:HN_ADNI / tab:HN_NACC)
%% ================================================================
function write_latex_table(metrics, fname, cohort, bw, est_tag, has_pcal, all_loess)
    c = char(cohort);
    if all_loess
        curves_txt = sprintf(['All curves, global and of each subgroup, are estimated with %s ' ...
            '($bw=%g$~years). '], est_tag, bw);
    else
        curves_txt = sprintf(['Global: polynomial trend. Amyloid and vascular: curve of each ' ...
            'subgroup (%s with $bw=%g$~years, or polynomial). '], est_tag, bw);
    end
    caption = sprintf(['Metrics of the trajectories of Fig.~\\ref{M-fig:HN_%s} of the main article, %s. ' ...
        '$t$: years relative to the transition age estimated by the DPM ($t=0$). ' ...
        '$t_{\\rm start}$, $t_{\\rm end}$: limits of the window where each curve is estimated. ' ...
        '$y_{\\rm start}$, $y_{\\rm end}$: value of the curve at those limits. ' ...
        'Slope: $(y_{\\rm end}-y_{\\rm start})/(t_{\\rm end}-t_{\\rm start})$, in units/year (MEM, EF), CL/year (Centiloids) or ml/year (WMH). ' ...
        '%s' ...
        '$p$: permutation test on the mean absolute difference between the two curves of each comparison.%s ' ...
        '$^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$.}'], c, c, curves_txt, ...
        ternary(has_pcal, ' $p_{\rm cal}$: simulation-calibrated $p$-value.', ''));
    caption = caption(1:end-1);   % the closing brace is added by the writer
    write_metrics_table_es(metrics, fname, caption, sprintf('tab:HN_%s', c), ...
        'Analysis', {'Global', 'Amyloid', 'Vascular'}, has_pcal);
end


function s = fmt_pval(p)
    if p < 0.001
        s = '$<$0.001$^{***}$';
    elseif p < 0.01
        s = sprintf('%.3f$^{**}$', p);
    elseif p < 0.05
        s = sprintf('%.3f$^{*}$', p);
    else
        s = sprintf('%.3f', p);
    end
end


function s = escape_tex(s)
    s = strrep(s, '_', '\_');
    s = strrep(s, 'A-', 'A$-$');
    s = strrep(s, 'V-', 'V$-$');
    s = strrep(s, 'A+', 'A$+$');
    s = strrep(s, 'V+', 'V$+$');
end


%% ================================================================
%%  GENERIC HELPERS (plotting / smoothing)
%% ================================================================

function tr = subgrp_range(data, range_p, t_range)
    valid = ~isnan(data.Time_to_Event);
    if sum(valid) < 5, tr = t_range; return; end
    u_ids  = unique(data.SubjID(valid));
    t_subj = arrayfun(@(id) median(data.Time_to_Event( ...
                 data.SubjID==id & valid)), u_ids);
    t_subj = t_subj(~isnan(t_subj));
    if numel(t_subj) < 3, tr = t_range; return; end
    t_lo = max(prctile(t_subj, range_p),       t_range(1));
    t_hi = min(prctile(t_subj, 100-range_p),   t_range(2));
    if t_hi <= t_lo, tr = t_range; return; end
    if (t_hi-t_lo) >= 0.90*(t_range(2)-t_range(1))
        tr = t_range;
    else
        tr = [t_lo, t_hi];
    end
end


function h = plot_lines_local_h(data, varName, color)
    h = gobjects(1);
    if isempty(data), return; end
    u_ids = unique(data.SubjID);
    for i = 1:length(u_ids)
        subj = sortrows(data(data.SubjID==u_ids(i),:),'Time_to_Event');
        mask = ~isnan(subj.(varName));
        if sum(mask) >= 2
            plot(subj.Time_to_Event(mask), subj.(varName)(mask), ...
                 'Color',[color 0.15],'LineWidth',0.5,'HandleVisibility','off');
        end
    end
    h = plot(nan, nan, '-', 'Color', color, 'LineWidth', 1.5, 'HandleVisibility','on');
end


% Point cloud with an optional connecting line, for sparse longitudinal
% variables (NACC CENTILOIDS / WMH). Copied verbatim from
% plot_naturalhistory_NACC_v2_1.m.
%   - 1 valid measurement  -> isolated marker, alpha=0.30
%   - >=2 measurements     -> connected markers, alpha=0.45
function h = plot_scatter_local_h(data, varName, color)
    h = gobjects(1);
    if isempty(data) || ~ismember(varName, data.Properties.VariableNames), return; end
    u_ids = unique(data.SubjID);
    for i = 1:length(u_ids)
        subj = sortrows(data(data.SubjID==u_ids(i),:),'Time_to_Event');
        mask = ~isnan(subj.(varName)) & ~isnan(subj.Time_to_Event);
        n_pts = sum(mask);
        if n_pts == 0, continue; end
        t_pts = subj.Time_to_Event(mask);
        y_pts = subj.(varName)(mask);
        if n_pts == 1
            scatter(t_pts, y_pts, 18, color, 'filled', ...
                    'MarkerFaceAlpha', 0.30, 'HandleVisibility','off');
        else
            plot(t_pts, y_pts, '-', 'Color', [color 0.25], 'LineWidth', 0.7, ...
                 'HandleVisibility','off');
            scatter(t_pts, y_pts, 18, color, 'filled', ...
                    'MarkerFaceAlpha', 0.45, 'HandleVisibility','off');
        end
    end
    h = scatter(nan, nan, 20, color, 'filled', 'HandleVisibility','on');
end


function [h, y_plot, x_plot] = plot_trend_poly_out(data, varName, color, width, ...
                                  poly_degree, density_trim_p, t_range, ...
                                  show_CI, n_boot, n_query)
    h = gobjects(1); y_plot = []; x_plot = [];
    if ~ismember(varName, data.Properties.VariableNames), return; end
    validMask = ~isnan(data.(varName)) & ~isnan(data.Time_to_Event);
    x_all = data.Time_to_Event(validMask);
    y_all = data.(varName)(validMask);
    if length(x_all) < 10, return; end
    t_lo = max(t_range(1), prctile(x_all, density_trim_p));
    t_hi = min(t_range(2), prctile(x_all, 100-density_trim_p));
    if t_hi <= t_lo, return; end
    mask = x_all >= t_lo & x_all <= t_hi;
    x = x_all(mask); y = y_all(mask);
    if length(x) < poly_degree+2, return; end
    x_plot = linspace(t_lo, t_hi, n_query);
    p      = polyfit(x, y, poly_degree);
    y_plot = polyval(p, x_plot);
    if show_CI
        y_res = y - polyval(p, x);
        Y_boot = zeros(n_boot, n_query);
        for b = 1:n_boot
            yb = polyval(p,x) + y_res(randi(length(y_res),length(y_res),1));
            Y_boot(b,:) = polyval(polyfit(x,yb,poly_degree), x_plot);
        end
        fill([x_plot,fliplr(x_plot)],[prctile(Y_boot,2.5,1),fliplr(prctile(Y_boot,97.5,1))], ...
             color,'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
    end
    h = plot(x_plot, y_plot, 'Color', color, 'LineWidth', width, 'HandleVisibility','on');
end


function [h, y_plot, x_plot] = plot_trend_loess_sw_out(data, varName, color, width, ...
                                  bw, t_range, show_CI, n_boot, n_query, smoother)
    h = gobjects(1); y_plot = []; x_plot = [];
    if ~ismember(varName, data.Properties.VariableNames), return; end
    validMask = ~isnan(data.(varName)) & ~isnan(data.Time_to_Event);
    if sum(validMask) < 10, return; end
    t_all  = data.Time_to_Event(validMask);
    y_all  = data.(varName)(validMask);
    id_all = data.SubjID(validMask);
    u_ids  = unique(id_all);
    n_subj = numel(u_ids);
    if n_subj < 5, return; end
    t_lo = max(t_range(1), prctile(t_all,  5));
    t_hi = min(t_range(2), prctile(t_all, 95));
    if t_hi <= t_lo, return; end
    x_plot = linspace(t_lo, t_hi, n_query);
    y_plot = smoother(t_all, y_all, id_all, u_ids, x_plot, bw);
    if all(isnan(y_plot)), y_plot = []; x_plot = []; return; end
    if show_CI
        Y_boot = nan(n_boot, n_query);
        for b = 1:n_boot
            idx_b = randi(n_subj, n_subj, 1);
            ids_b = u_ids(idx_b);
            mask_b = ismember(id_all, ids_b);
            Y_boot(b,:) = smoother(t_all(mask_b), y_all(mask_b), ...
                                   id_all(mask_b), ids_b, x_plot, bw);
        end
        fill([x_plot,fliplr(x_plot)], ...
             [prctile(Y_boot,2.5,1), fliplr(prctile(Y_boot,97.5,1))], ...
             color,'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
    end
    h = plot(x_plot, y_plot, 'Color', color, 'LineWidth', width, 'HandleVisibility','on');
end


function y_hat = sw_loess_eval(t, y, ids, u_ids, x_plot, bw)
    % Subject-weighted Nadaraya-Watson smoother (sw-LOESS), exactly as in
    % Eqs. (2)-(4) of the manuscript (Sec. 3.3.1):
    %   k_ij(t)   = exp(-(t - t*_ij)^2 / (2 bw^2))               (2)
    %   ybar_i(t) = sum_j k_ij(t) y_ij / sum_j k_ij(t)           (3)
    %   yhat(t)   = sum_i w_i(t) ybar_i(t) / sum_i w_i(t),
    %               w_i(t) = max_j k_ij(t)                       (4)
    % The subject weight w_i(t) depends only on the subject's visit
    % closest to t, not on how many visits it has. The previous version
    % used w_i(t) = sum_j k_ij(t); with that choice num/denom collapses
    % to sum_ij k_ij y_ij / sum_ij k_ij, i.e. the classical
    % observation-level Nadaraya-Watson estimator, not a subject-weighted
    % one.
    n_q = numel(x_plot);
    y_hat = nan(1, n_q);
    n_s = numel(u_ids);
    t_cell = cell(n_s,1); y_cell = cell(n_s,1);
    for si = 1:n_s
        m = ids == u_ids(si);
        t_cell{si} = t(m); y_cell{si} = y(m);
    end
    for qi = 1:n_q
        t0 = x_plot(qi); num = 0; denom = 0;
        for si = 1:n_s
            k_s = exp(-0.5*((t_cell{si}-t0)/bw).^2);   % k_ij(t), Eq. (2)
            if isempty(k_s), continue; end
            w_s = max(k_s);                              % w_i(t) = max_j k_ij(t), Eq. (4)
            if w_s < 1e-12, continue; end
            ybar_s = sum(k_s.*y_cell{si}) / sum(k_s);    % ybar_i(t), Eq. (3)
            num   = num   + w_s*ybar_s;
            denom = denom + w_s;
        end
        if denom > 1e-12, y_hat(qi) = num/denom; end
    end
end


function m = make_metric(g, row, col, var, subgrp, xp, yp, pval, n_subj, pcal)
    m.group   = g;   m.row = row;  m.col = col;
    m.var     = var; m.subgrp = subgrp;
    m.t_start = xp(1);   m.t_end = xp(end);
    m.y_start = yp(1);   m.y_end = yp(end);
    m.slope   = (yp(end)-yp(1))/(xp(end)-xp(1));
    m.pval    = pval;    m.n_subj = n_subj;    m.pcal = pcal;
end


function build_legend(handles, labels)
    % See the header note: no get(h,'Type') call here, unlike the
    % version originally in plot_tau_dif_batch_v2k.m. It errors on a
    % GraphicsPlaceholder, and isvalid(h) does not reliably prevent
    % reaching it first when a subgroup is empty.
    valid = arrayfun(@(h) ~isempty(h) && isvalid(h) && ...
                         ~isa(h,'matlab.graphics.GraphicsPlaceholder'), handles);
    h_ok = handles(valid); l_ok = labels(valid);
    if isempty(h_ok), return; end
    legend(h_ok, l_ok, 'Location','best','FontSize',7);
end


%% ================================================================
%%  v2 HELPERS: options, p-value tables, figure saving, smoothers
%% ================================================================

function opt = parse_options(varargin)
    opt = struct('estimator', 'linear', 'pFile', '', 'pcalFile', '', ...
                 'saveDir', '', 'figVersion', 'v3', 'allLoess', true, ...
                 'dataDir', 'data', 'tauDir', fullfile('data', 'tau'), ...
                 'wmhThreshold', 10, 'outDir', '.');
    if mod(numel(varargin), 2) ~= 0
        error('parse_options:pairs', 'Options must be name-value pairs.');
    end
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~isfield(opt, name)
            error('parse_options:unknown', 'Unknown option ''%s''.', name);
        end
        opt.(name) = varargin{k+1};
    end
    opt.estimator = lower(char(opt.estimator));
    if ~ismember(opt.estimator, {'linear', 'const'})
        error('parse_options:estimator', 'estimator must be ''linear'' or ''const''.');
    end
    opt.allLoess = logical(opt.allLoess);
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end

% Permutation p-values of bw_sensitivity (same permutations as Table 4
% of paper A). Only the rows of the bandwidth in use are kept.
function PT = read_p_table(file, estimator, bw)
    PT = [];
    if isempty(file), return; end
    [~, nm] = fileparts(file);
    isll = numel(nm) > 3 && strcmp(nm(end-2:end), '_ll');
    if isll ~= strcmp(estimator, 'linear')
        warning('read_p_table:estimator', ...
            '%s does not match estimator ''%s'' (''_ll'' = linear).', file, estimator);
    end
    if ~exist(file, 'file')
        error('read_p_table:file', 'pFile not found: %s (current folder: %s).', file, pwd);
    end
    T = readtable(file, 'TextType', 'string');
    PT = T(abs(T.bw - bw) < 1e-9, {'transition', 'analysis', 'var', 'p'});
end

% Calibrated p-values of sim_framework (sim_results.csv, design
% 'calibrated'; scenario = 'COHORT|TRANS|ANALYSIS|VAR').
function PCT = read_pcal_table(file, estimator, bw)
    PCT = [];
    if isempty(file), return; end
    if ~exist(file, 'file')
        error('read_pcal_table:file', 'pcalFile not found: %s (current folder: %s).', file, pwd);
    end
    T = readtable(file, 'TextType', 'string');
    est = ternary(strcmp(estimator, 'linear'), sprintf('swll_bw%g', bw), sprintf('sw_bw%g', bw));
    PCT = T(T.design == "calibrated" & T.estimator == est, {'scenario', 'pCal'});
end

function p = lookup_p(PT, tcode, analysis, vn)
    p = NaN;
    if isempty(PT), return; end
    m = PT.transition == string(tcode) & PT.analysis == string(analysis) & PT.var == string(vn);
    if any(m), p = PT.p(find(m, 1)); end
end

function p = lookup_pcal(PCT, cohort, tcode, analysis, vn)
    p = NaN;
    if isempty(PCT), return; end
    key = sprintf('%s|%s|%s|%s', char(cohort), char(tcode), char(analysis), char(vn));
    m = PCT.scenario == string(key);
    if any(m), p = PCT.pCal(find(m, 1)); end
end

function save_figure(fh, file)
    d = fileparts(file);
    if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end
    if exist('exportgraphics', 'file') == 2
        exportgraphics(fh, file, 'Resolution', 300);
    else
        print(fh, file, '-dpng', '-r300');
    end
    fprintf('  Figure saved to %s\n', file);
end

function s = fmt_p_es(p)
    % p-value with stars, as in tab:HN_ADNI / tab:HN_NACC of the supplement
    if isnan(p)
        s = 'n.e.';
    elseif p < 0.001
        s = '$<$0.001$^{***}$';
    elseif p < 0.01
        s = sprintf('%.3f$^{**}$', p);
    elseif p < 0.05
        s = sprintf('%.3f$^{*}$', p);
    else
        s = sprintf('%.3f', p);
    end
end

function s = fmt_pcal_es(p)
    if isnan(p)
        s = '---';
    elseif p < 0.0021
        s = '$\leq$0.002';
    else
        s = sprintf('%.3f', p);
    end
end

function s = var_es(vn)
    s = strrep(char(vn), 'CENTILOIDS', 'Centiloids');
end

% LaTeX table, one block per comparison (group x row x column),
% p-value (and p_cal) on the first line of the block via \multirow.
function write_metrics_table_es(metrics, fname, caption, label, col2_name, row_names, has_pcal)
    fid = fopen(fname, 'w');
    ncol = 10 + has_pcal;
    fprintf(fid, '%% Auto-generated (%s)\n\\begin{table}[htbp]\n\t\\centering\n', datestr(now, 'yyyy-mm-dd'));
    fprintf(fid, '\t\\caption{%s}\n\t\\label{%s}\n', caption, label);
    fprintf(fid, '\t\\setlength{\\tabcolsep}{4pt}\n\t\\renewcommand{\\arraystretch}{1.15}\n\t\\tiny\n');
    fprintf(fid, '\t\\begin{tabular}{llllrrrrr%s}\n\t\t\\toprule\n', ternary(has_pcal, 'll', 'l'));
    fprintf(fid, ['\t\tTransition & %s & Variable & Subgroup & $t_{\\rm start}$ & $y_{\\rm start}$\n' ...
                  '\t\t& $t_{\\rm end}$ & $y_{\\rm end}$ & Slope & $p$%s \\\\\n\t\t\\midrule\n'], ...
            col2_name, ternary(has_pcal, ' & $p_{\rm cal}$', ''));
    trans_tex = {'CU$\to$MCI', 'MCI$\to$Dem.'};
    ord = zeros(numel(metrics), 4);
    for k = 1:numel(metrics)
        sg = metrics(k).subgrp;
        ord(k,:) = [metrics(k).group, metrics(k).row, metrics(k).col, ...
                    any(strcmp(sg, {'A+', 'V+'}))];
    end
    [~, sidx] = sortrows(ord);  metrics = metrics(sidx);  ord = ord(sidx, :);
    blocks = unique(ord(:, 1:3), 'rows', 'stable');
    for b = 1:size(blocks, 1)
        idx = find(ismember(ord(:, 1:3), blocks(b, :), 'rows'))';
        if b > 1
            if blocks(b, 1) ~= blocks(b-1, 1)
                fprintf(fid, '\t\t\\midrule[0.8pt]\n');
            elseif blocks(b, 2) > 0 && ~(blocks(b, 2) == 1 && blocks(b-1, 2) == 1)
                fprintf(fid, '\t\t\\cmidrule{3-%d}\n', ncol);
            end
        end
        for j = 1:numel(idx)
            m = metrics(idx(j));
            if isnan(m.slope), sl = '---'; else, sl = sprintf('%.4f', m.slope); end
            if strcmp(m.subgrp, 'Global')
                ps = '---';  pcs = '---';
            elseif j == 1
                ps  = sprintf('\\multirow{%d}{*}{%s}', numel(idx), fmt_p_es(m.pval));
                pcs = sprintf('\\multirow{%d}{*}{%s}', numel(idx), fmt_pcal_es(m.pcal));
            else
                ps = '';  pcs = '';
            end
            sg = strrep(strrep(strrep(strrep(m.subgrp, 'A-', 'A$-$'), 'A+', 'A$+$'), 'V-', 'V$-$'), 'V+', 'V$+$');
            fprintf(fid, '\t\t%s & %s & %s & %s & %.2f & %.3f & %.2f & %.3f & %s & %s', ...
                trans_tex{m.group}, row_names{m.row}, var_es(m.var), sg, ...
                m.t_start, m.y_start, m.t_end, m.y_end, sl, ps);
            if has_pcal, fprintf(fid, ' & %s', pcs); end
            fprintf(fid, ' \\\\\n');
        end
    end
    fprintf(fid, '\t\t\\bottomrule\n\t\\end{tabular}\n\\end{table}\n');
    fclose(fid);
end

function y_hat = sw_loess_ll_eval(t, y, ids, u_ids, x_plot, bw)
    % Subject-weighted LOCAL-LINEAR smoother (sw-LOESS, Eqs. 5-6 of paper
    % A, v3). Same subject weight as sw_loess_eval, w_i(t) = max_j k_ij(t),
    % shared among the subject's visits, a_ij(t) = w_i(t) k_ij(t) /
    % sum_j k_ij(t). At each t a weighted line beta0 + beta1 (t*_ij - t)
    % is fitted, and the curve is beta0. With beta1 = 0 it reduces to
    % sw_loess_eval. Duplicated ids in u_ids (subject bootstrap) count
    % that subject as many times as it appears.
    n_q = numel(x_plot);
    y_hat = nan(1, n_q);
    S0 = zeros(1, n_q); S1 = S0; S2 = S0; T0 = S0; T1 = S0;
    x_plot = x_plot(:)';
    for si = 1:numel(u_ids)
        m = ids == u_ids(si);
        if ~any(m), continue; end
        dt = t(m) - x_plot;                         % visits x n_q
        K  = exp(-0.5 * (dt / bw).^2);
        w  = max(K, [], 1);  w(w < 1e-12) = 0;
        a  = K .* (w ./ max(sum(K, 1), realmin));
        ys = y(m);
        S0 = S0 + sum(a, 1);           S1 = S1 + sum(a .* dt, 1);
        S2 = S2 + sum(a .* dt.^2, 1);  T0 = T0 + sum(a .* ys, 1);
        T1 = T1 + sum(a .* dt .* ys, 1);
    end
    den = S0 .* S2 - S1.^2;
    ok = den > 1e-10 * max(S0.^2, realmin);
    y_hat(ok) = (S2(ok) .* T0(ok) - S1(ok) .* T1(ok)) ./ den(ok);
end
