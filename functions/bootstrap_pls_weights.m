function boot_out = bootstrap_pls_weights(X, Y, observed_w, B, store_full, verbose_interval)
% bootstrap_pls_weights
%   X: (nSubj x nEdges) predictor matrix (already zscored)
%   Y: (nSubj x 1) behavioral vector (already mean-imputed; raw or z)
%   observed_w: (nEdges x 1) weights from original (full-sample) PLS
%   B: number of bootstrap resamples (e.g., 1000)
%   store_full: logical; if true, store full boot weight matrix
%   verbose_interval: print progress every this many iterations
%
% Returns struct with fields described in main text.

    if nargin < 5 || isempty(store_full), store_full = false; end
    if nargin < 6 || isempty(verbose_interval), verbose_interval = max(1,floor(B/10)); end

    [nSubj, nEdges] = size(X);

    % Online accumulators (Welford)
    w_mean = zeros(nEdges,1);
    w_M2   = zeros(nEdges,1);
    sign_same = zeros(nEdges,1); % count replicates having same sign as observed
    obs_sign = sign(observed_w);
    obs_sign(obs_sign==0) = 1;   % treat zeros as positive for stability counting

    if store_full
        % This can be large: nEdges x B (≈ 4005 x 1000 ≈ 32 MB double)
        W_all = zeros(nEdges, B, 'single'); % single precision to save memory
    end

    for b = 1:B
        idx = randsample(nSubj, nSubj, true);
        Xb  = X(idx,:);
        Yb  = Y(idx);
        % Fit 1-component PLS
        try
            [~,~,~,~,~,~,~,stats_b] = plsregress(Xb, zscore(Yb), 1);
            w_b = stats_b.W;  % nEdges x 1
        catch
            % Fallback: if numerical issue, treat as zeros
            w_b = zeros(nEdges,1);
        end

        % Online mean / variance
        delta = w_b - w_mean;
        w_mean = w_mean + delta / b;
        w_M2   = w_M2 + delta .* (w_b - w_mean);

        % Sign stability
        sign_same = sign_same + (sign(w_b) == obs_sign);

        % Store full matrix if requested
        if store_full
            W_all(:, b) = single(w_b);
        end

        if mod(b, verbose_interval) == 0
            fprintf('  Bootstrap %d/%d completed.\n', b, B);
        end
    end

    SE = sqrt(w_M2 / max(1, B - 1));
    SE(SE==0) = Inf; % avoid division by zero; ratio -> 0

    BR = observed_w ./ SE; % bootstrap ratio
    p_norm = 2*normcdf(-abs(BR), 0, 1);

    % Empirical sign proportion
    signProp = sign_same / B;

    % Flag stable weights
    flag_BR3 = abs(BR) > 3;

    boot_out = struct();
    boot_out.boot_n        = B;
    boot_out.boot_mean     = w_mean;
    boot_out.boot_SE       = SE;
    boot_out.boot_BR       = BR;
    boot_out.boot_signProp = signProp;
    boot_out.boot_flag_BR3 = flag_BR3;
    boot_out.boot_p_norm   = p_norm;
    % FDR across edges
    [~, sortIdx] = sort(p_norm);
    ranks = zeros(size(sortIdx));
    ranks(sortIdx) = 1:numel(p_norm);
    boot_out.boot_p_FDR = p_norm .* numel(p_norm) ./ ranks;
    boot_out.boot_p_FDR(boot_out.boot_p_FDR > 1) = 1;

    if store_full
        % 95% percentile CI
        boot_out.boot_CI95_lo = quantile(double(W_all), 0.025, 2);
        boot_out.boot_CI95_hi = quantile(double(W_all), 0.975, 2);
        % Optionally discard W_all to save disk; or keep:
        % boot_out.boot_weights_full = W_all;  % COMMENT OUT if too large
    end
end