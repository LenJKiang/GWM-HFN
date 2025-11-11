% Batch_Calculate_ICC_and_Compare.m
% This script calculates the ICC (intraclass correlation coefficient) for GM-GM, GM-WM, and GWM-HFN networks
% for both SLIM and BNU3 datasets. It also compares the ICC results of GM-GM and GWM-HFN using paired t-tests
% and bootstrapping, following the conventions used in Seguin 2022 and previous analyses.
% Author: Zeqiang Linli, 2025-08-03

clear; clc;

%% ====== Dataset and Path Configuration ======
dbs = {'SLIM','BNU3'};
db_dir.SLIM  = 'E:\Neuroimage\MyProject\GMWM_Network\Data\SLIM\';
db_dir.BNU3  = 'E:\Neuroimage\MyProject\GMWM_Network\Data\BNU3\';
db_sessions.SLIM = 2;   % SLIM: 2 sessions
db_sessions.BNU3 = 3;   % BNU3: 3 sessions

% Use internal field names expected downstream (no underscore)
network_types = {'GMGM','GWMHFN'};

ICC_results = struct();

for d = 1:numel(dbs)
    db = dbs{d};
    n_session = db_sessions.(db);
    base_dir = db_dir.(db);

    fprintf('\n=== Calculating ICC for %s dataset ===\n', db);

    mats = struct();
    ids_all = struct();

    % ---------- Minimal, robust data loader (SLIM vs BNU3) ----------
    for t = 1:n_session
        for n = 1:numel(network_types)
            net = network_types{n};          % internal name, e.g. 'GWMHFN' or 'GMGM'
            % file_base is the token used in filenames/variable names on disk
            file_base = net;
            if strcmp(net, 'GWMHFN')
                file_base = 'GWM_HFN';      % disk uses underscore for this network
            end

            if strcmpi(db, 'SLIM')
                % Prefer the "common" file (paired common subjects). Fallback to "all" file.
                fn_common = fullfile(base_dir, sprintf('%s_time%d_common.mat', file_base, t));
                fn_all    = fullfile(base_dir, sprintf('%s_time%d.mat',      file_base, t));
                if exist(fn_common, 'file')
                    tmp = load(fn_common);
                    pref_var = sprintf('time%d_%s_common', t, strrep(file_base,'_',''));
                elseif exist(fn_all, 'file')
                    tmp = load(fn_all);
                    pref_var = sprintf('time%d_%s', t, strrep(file_base,'_',''));
                else
                    error('SLIM expected file not found for %s time%d: tried %s and %s', file_base, t, fn_common, fn_all);
                end
            else
                % BNU3: expect time{t}_ files (all subjects)
                fn = fullfile(base_dir, sprintf('%s_time%d.mat', file_base, t));
                if ~exist(fn, 'file')
                    error('BNU3 expected file not found: %s', fn);
                end
                tmp = load(fn);
                pref_var = sprintf('time%d_%s', t, strrep(file_base,'_',''));
            end

            % Assign: if preferred variable exists in the .mat, use it; otherwise pick first numeric matrix-like var
            if isfield(tmp, pref_var)
                mats.(sprintf('time%d_%s', t, net)) = tmp.(pref_var);
            else
                % fallback: take first numeric variable with ndims >= 2
                fnames = fieldnames(tmp);
                found = false;
                for fi = 1:numel(fnames)
                    cand = tmp.(fnames{fi});
                    if isnumeric(cand) && ndims(cand) >= 2
                        mats.(sprintf('time%d_%s', t, net)) = cand;
                        found = true;
                        break;
                    end
                end
                if ~found
                    error('No suitable numeric variable found in file for db=%s, time=%d, net=%s', db, t, net);
                end
            end
        end

        % Load subject IDs (saved as time{t}_subject_ids.mat with variable time{t}_subject_ids,
        % but fallback to 'ids' or first variable)
        id_file = fullfile(base_dir, sprintf('time%d_subject_ids.mat', t));
        if exist(id_file, 'file')
            idS = load(id_file);
            preferred_id_var = sprintf('time%d_subject_ids', t);
            if isfield(idS, preferred_id_var)
                ids_all.(sprintf('time%d', t)) = idS.(preferred_id_var);
            elseif isfield(idS, 'ids')
                ids_all.(sprintf('time%d', t)) = idS.ids;
            else
                fn = fieldnames(idS);
                if ~isempty(fn)
                    ids_all.(sprintf('time%d', t)) = idS.(fn{1});
                else
                    ids_all.(sprintf('time%d', t)) = [];
                end
            end
        else
            ids_all.(sprintf('time%d', t)) = [];
        end
    end
    % ---------- end loader ----------

    % Note: SLIM saved "common" versions already (if available), so no additional alignment is needed.

    %% --- Calculate ICC for each network type ---
    for n = 1:numel(network_types)
        net = network_types{n};
        fprintf('\n--- %s: %s ---\n', db, net);

        if n_session == 2
            % two timepoints
            data1 = mats.(sprintf('time1_%s', net)); % [node x node x subj]
            data2 = mats.(sprintf('time2_%s', net));
            ICC_mat = calculate_icc_matrix(data1, data2);
        elseif n_session == 3
            % three timepoints
            data1 = mats.(sprintf('time1_%s', net));
            data2 = mats.(sprintf('time2_%s', net));
            data3 = mats.(sprintf('time3_%s', net));
            ICC_mat = calculate_icc_matrix_three_sessions(data1, data2, data3);
        end

        % Extract lower triangle values (excluding diagonal) for summary
        lower_idx = find(tril(ones(size(ICC_mat)), -1));
        ICC_vec = ICC_mat(lower_idx);

        ICC_results.(db).(net).ICC_matrix = ICC_mat;
        ICC_results.(db).(net).ICC_vec = ICC_vec;
        ICC_results.(db).(net).mean_ICC = mean(ICC_vec, 'omitnan');
        ICC_results.(db).(net).std_ICC = std(ICC_vec, 'omitnan');
    end

    %% --- Compare GM-GM and GWM-HFN ICC distributions ---
    % Note: internal field name is 'GWMHFN' (no underscore)
    ICC1 = ICC_results.(db).GMGM.ICC_vec;
    ICC2 = ICC_results.(db).GWMHFN.ICC_vec;
    [~, p_ttest, ~, stats] = ttest(ICC1, ICC2);

    % Bootstrap confidence interval for mean difference
    n_boot = 10000;
    boot_diff = zeros(n_boot,1);
    n_edge = length(ICC1);
    rng(12345); % for reproducibility
    for b = 1:n_boot
        idx = randsample(n_edge, n_edge, true);
        boot_diff(b) = mean(ICC1(idx) - ICC2(idx));
    end
    ci = prctile(boot_diff, [2.5 97.5]);
    p_boot = 2 * min(mean(boot_diff <= 0), mean(boot_diff > 0));

    ICC_results.(db).Comparison.GMGM_vs_GWMHFN.ttest_p = p_ttest;
    ICC_results.(db).Comparison.GMGM_vs_GWMHFN.ttest_t = stats.tstat;
    ICC_results.(db).Comparison.GMGM_vs_GWMHFN.mean_diff = mean(ICC1 - ICC2);
    ICC_results.(db).Comparison.GMGM_vs_GWMHFN.ci_95 = ci;
    ICC_results.(db).Comparison.GMGM_vs_GWMHFN.p_boot = p_boot;
    ICC_results.(db).Comparison.GMGM_vs_GWMHFN.boot_diff = boot_diff;  % save bootstrap vector

    fprintf('\n--- Edgewise ICC comparison: GM-GM vs GWM-HFN (%s) ---\n', db);
    fprintf('Mean ICC (GM-GM): %.4f, Mean ICC (GWM-HFN): %.4f\n', mean(ICC1), mean(ICC2));
    fprintf('Paired t-test: t = %.2f, p = %.4g\n', stats.tstat, p_ttest);
    fprintf('Mean difference: %.4f, 95%% CI = [%.4f, %.4f], bootstrap p = %.4g\n', ...
        mean(ICC1 - ICC2), ci(1), ci(2), p_boot);
end

%% ====== Display and Save Results ======
fprintf('\n======= ICC Summary =======\n');
for d = 1:numel(dbs)
    db = dbs{d};
    fprintf('\n=== %s ===\n', db);
    for n = 1:numel(network_types)
        net = network_types{n};
        r = ICC_results.(db).(net);
        fprintf('%-8s: Mean ICC = %.4f, SD = %.4f\n', net, r.mean_ICC, r.std_ICC);
    end
    % Print comparison
    cmp = ICC_results.(db).Comparison.GMGM_vs_GWMHFN;
    fprintf('GMGM vs GWMHFN: t = %.2f, p = %.4g, diff = %.4f, 95%% CI = [%.4f, %.4f], boot p = %.4g\n', ...
        cmp.ttest_t, cmp.ttest_p, cmp.mean_diff, cmp.ci_95(1), cmp.ci_95(2), cmp.p_boot);
end

save('ICC_Results_All.mat', 'ICC_results', '-v6');
fprintf('\nAll ICC values computed and saved to ICC_Results_All.mat\n');