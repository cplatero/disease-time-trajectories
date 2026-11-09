function y_hat = sw_loess_ll(t, y, ids, u_ids, x_plot, bw)
    % SW_LOESS_LL  Subject-weighted local-linear smoother (sw-LOESS).
    %   y_hat = sw_loess_ll(t, y, ids, u_ids, x_plot, bw)
    %   t, y, ids: visits (disease time, marker, subject); u_ids: subjects
    %   used; x_plot: query points; bw: bandwidth (years).
    %   Public copy of sw_loess_ll_eval (plot_kinematics.m), for the tutorial.
    % (Eqs. 5-6 of paper
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

% Authorship: Carlos Platero (UPM), with the collaboration of Claude
% (Anthropic).
