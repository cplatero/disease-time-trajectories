function tr = coverage_range(data, range_p, t_range)
    % COVERAGE_RANGE  Coverage range of a subgroup on the disease-time axis:
    %   percentiles range_p and 100-range_p of the subject medians, within
    %   t_range. Public copy of subgrp_range (plot_kinematics.m), for the tutorial.
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

% Authorship: Carlos Platero (UPM), with the collaboration of Claude
% (Anthropic).
