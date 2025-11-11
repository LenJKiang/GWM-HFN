function q = mafdr(p, method, bh_correction)
    if nargin < 2, method = 'pdep'; end
    if nargin < 3, bh_correction = false; end
    p = p(:);
    valid_idx = ~isnan(p);
    p_valid = p(valid_idx);
    q = NaN(size(p));
    if isempty(p_valid), q(valid_idx) = p_valid; return; end
    m = length(p_valid);
    [sorted_p, sort_idx] = sort(p_valid);
    inv_idx(sort_idx) = 1:m;
    if strcmpi(method, 'pdep') || strcmpi(method, 'bhfdr') || bh_correction
        q_vals = m * sorted_p ./ (1:m)';
    else
        cm = sum(1./(1:m));
        q_vals = m * cm * sorted_p ./ (1:m)';
    end
    q_vals = cummin(q_vals(end:-1:1));
    q_vals = q_vals(end:-1:1);
    q_vals = min(q_vals, 1);
    q(valid_idx) = q_vals(inv_idx);
end
