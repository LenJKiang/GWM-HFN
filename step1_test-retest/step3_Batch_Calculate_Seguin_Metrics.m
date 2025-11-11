% Batch_Calculate_Seguin_Metrics.m
% Batch calculation of Seguin 2022 metrics (Reliability / Uniformity / Identifiability)
% Minimal and robust data loader for SLIM / BNU3 naming conventions
% Author: adapted for Zeqiang Linli, 2025-10-29

clear; clc;

%% ====== Dataset and Path Configuration ======
dbs = {'SLIM','BNU3'};
db_dir.SLIM  = 'E:\Neuroimage\MyProject\GMWM_Network\Data\SLIM\';
db_dir.BNU3  = 'E:\Neuroimage\MyProject\GMWM_Network\Data\BNU3\';
db_sessions.SLIM = 2; % SLIM dataset has 2 sessions
db_sessions.BNU3 = 3; % BNU3 dataset has 3 sessions
network_types = {'GMGM','GWMHFN'}; % internal names used in downstream functions

results = struct();

for d = 1:numel(dbs)
    db = dbs{d};
    n_session = db_sessions.(db);
    base_dir = db_dir.(db);

    fprintf('\n=== Processing %s dataset (%d sessions) ===\n', db, n_session);

    mats = struct();
    ids_all = struct();

    % ----------------- Robust loader (minimal branching) -----------------
    for t = 1:n_session
        for n = 1:numel(network_types)
            net = network_types{n};
            file_base = net;
            if strcmp(net, 'GWMHFN'); file_base = 'GWM_HFN'; end

            if strcmpi(db, 'SLIM')
                % SLIM: prefer common, then all, then plain
                fn_common = fullfile(base_dir, sprintf('%s_time%d_common.mat', file_base, t));
                fn_all    = fullfile(base_dir, sprintf('%s_time%d_all.mat',    file_base, t));
                fn_plain  = fullfile(base_dir, sprintf('%s_time%d.mat',        file_base, t));
                if exist(fn_common,'file')
                    tmp = load(fn_common);
                    cand_vars = {sprintf('time%d_%s_common', t, strrep(file_base,'_','')), ...
                        sprintf('%s_time%d_common', file_base, t)};
                    mats.(sprintf('time%d_%s', t, net)) = pick_var_or_fallback(tmp, cand_vars);
                elseif exist(fn_all,'file')
                    tmp = load(fn_all);
                    cand_vars = {sprintf('time%d_%s', t, strrep(file_base,'_','')), ...
                        sprintf('%s_time%d_all', file_base, t), ...
                        sprintf('%s_time%d', file_base, t)};
                    mats.(sprintf('time%d_%s', t, net)) = pick_var_or_fallback(tmp, cand_vars);
                elseif exist(fn_plain,'file')
                    tmp = load(fn_plain);
                    cand_vars = {sprintf('time%d_%s', t, strrep(file_base,'_','')), sprintf('%s_time%d', file_base, t)};
                    mats.(sprintf('time%d_%s', t, net)) = pick_var_or_fallback(tmp, cand_vars);
                else
                    error('SLIM: Missing file for %s time%d (tried common/all/plain)', file_base, t);
                end
            else
                % BNU3: expect file_base_time{t}.mat with variable time{t}_{file_base}
                fn = fullfile(base_dir, sprintf('%s_time%d.mat', file_base, t));
                if ~exist(fn,'file')
                    error('BNU3: expected file not found: %s', fn);
                end
                tmp = load(fn);
                cand_vars = {sprintf('time%d_%s', t, file_base), sprintf('time%d_%s', t, strrep(file_base,'_','')), ...
                    sprintf('%s_time%d', file_base, t)};
                mats.(sprintf('time%d_%s', t, net)) = pick_var_or_fallback(tmp, cand_vars);
            end
        end

        % Load subject IDs
        id_file = fullfile(base_dir, sprintf('time%d_subject_ids.mat', t));
        if exist(id_file, 'file')
            idS = load(id_file);
            pref = sprintf('time%d_subject_ids', t);
            if isfield(idS, pref)
                ids_all.(sprintf('time%d', t)) = idS.(pref);
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
    % ----------------- end loader -----------------


    %% --- Calculate Seguin metrics for each network type ---
    for n = 1:numel(network_types)
        net = network_types{n};
        fprintf('\n--- %s: %s ---\n', db, net);

        if n_session == 2
            % two timepoints
            data1 = mats.(sprintf('time1_%s', net));
            data2 = mats.(sprintf('time2_%s', net));
            sim_mat = calculate_similarity_matrix_seguin(data1, data2);
            results.(db).(net) = calculate_seguin_metrics(sim_mat, net, db);
        elseif n_session == 3
            % three timepoints
            data1 = mats.(sprintf('time1_%s', net));
            data2 = mats.(sprintf('time2_%s', net));
            data3 = mats.(sprintf('time3_%s', net));
            results.(db).(net) = calculate_seguin_metrics_three_sessions(data1, data2, data3, net);
        end
    end
end

%% ======  Statistical comparison  ======
rng(12345); % for reproducibility
n_boot = 10000;

for db = {'SLIM','BNU3'}
    dbn = db{1};
    fprintf('\n=== Statistical comparison: %s GM-GM vs GWM-HFN ===\n', dbn);

    % 1. Paired t-test for reliability
    rel1 = results.(dbn).GMGM.intra_similarities;
    rel2 = results.(dbn).GWMHFN.intra_similarities;
    [~, p_rel, ~, stats_rel] = ttest(rel1, rel2);

    % 2. t-test for uniformity (usually unpaired, but can use ttest2)
    uni1 = results.(dbn).GMGM.inter_similarities;
    uni2 = results.(dbn).GWMHFN.inter_similarities;
    [~, p_uni, ~, stats_uni] = ttest2(uni1, uni2);

    % 3. Identifiability: report difference and bootstrap
    id1 = results.(dbn).GMGM.identifiability;
    id2 = results.(dbn).GWMHFN.identifiability;
    diff_id = id1 - id2;

    % Bootstrap (also save distribution)
    nsubj = length(rel1);

    boot_id1 = nan(n_boot,1);
    boot_id2 = nan(n_boot,1);
    boot_diff = nan(n_boot,1);

    for b = 1:n_boot
        idx = randsample(nsubj, nsubj, true);       
        rel1b = rel1(idx);
        rel2b = rel2(idx);

        uni1b = uni1(randsample(numel(uni1), numel(rel1b), true));
        uni2b = uni2(randsample(numel(uni2), numel(rel2b), true));

        % pooled std for network1 (rel1b vs uni1b)
        s_intra1 = std(rel1b, 0, 1);
        s_inter1 = std(uni1b, 0, 1);
        n_intra1 = numel(rel1b);
        n_inter1 = numel(uni1b);
        denom1 = ( (n_intra1-1)*s_intra1^2 + (n_inter1-1)*s_inter1^2 ) / (n_intra1 + n_inter1 - 2);
        s_pooled1 = sqrt(max(denom1, eps)); % 避免0

        % pooled std for network2 (rel2b vs uni2b)
        s_intra2 = std(rel2b, 0, 1);
        s_inter2 = std(uni2b, 0, 1);
        n_intra2 = numel(rel2b);
        n_inter2 = numel(uni2b);
        denom2 = ( (n_intra2-1)*s_intra2^2 + (n_inter2-1)*s_inter2^2 ) / (n_intra2 + n_inter2 - 2);
        s_pooled2 = sqrt(max(denom2, eps));

        % compute identifiability per network
        id1b = (mean(rel1b) - mean(uni1b)) / s_pooled1;
        id2b = (mean(rel2b) - mean(uni2b)) / s_pooled2;

        boot_id1(b) = id1b;
        boot_id2(b) = id2b;
        boot_diff(b) = id1b - id2b;
    end

    ci = prctile(boot_diff, [2.5 97.5]);
    p_boot = mean(boot_diff <= 0);

    results.(dbn).Comparison.GMGM_vs_GWMHFN.boot_diff = boot_diff;
    results.(dbn).Comparison.GMGM_vs_GWMHFN.boot_id1 = boot_id1;
    results.(dbn).Comparison.GMGM_vs_GWMHFN.boot_id2 = boot_id2;
    results.(dbn).Comparison.GMGM_vs_GWMHFN.boot_ci = ci;
    results.(dbn).Comparison.GMGM_vs_GWMHFN.p_boot = p_boot;
    results.(dbn).Comparison.GMGM_vs_GWMHFN.id_diff_point_est = diff_id;
    results.(dbn).Comparison.GMGM_vs_GWMHFN.id1_point_est = id1;
    results.(dbn).Comparison.GMGM_vs_GWMHFN.id2_point_est = id2;

    % Print results
    fprintf('Reliability: mean(GMGM)=%.4f, mean(GWM-HFN)=%.4f, t=%.2f, p=%.4g\n', ...
        mean(rel1), mean(rel2), stats_rel.tstat, p_rel);
    fprintf('Uniformity: mean(GMGM)=%.4f, mean(GWM-HFN)=%.4f, t=%.2f, p=%.4g\n', ...
        mean(uni1), mean(uni2), stats_uni.tstat, p_uni);
    fprintf('Identifiability: GMGM=%.4f, GWM-HFN=%.4f, diff=%.4f, 95%% CI=[%.4f, %.4f], bootstrap p=%.4g\n', ...
        id1, id2, diff_id, ci(1), ci(2), p_boot);
end

%% ====== Display and Save Results ======
fprintf('\n======= Seguin Reliability Metrics Summary =======\n');
for d = 1:numel(dbs)
    db = dbs{d};
    fprintf('\n=== %s ===\n', db);
    for n = 1:numel(network_types)
        net = network_types{n};
        r = results.(db).(net);
        fprintf('%-8s: Reliability = %.4f, Uniformity = %.4f, Identifiability = %.4f\n', ...
            net, r.reliability, r.uniformity, r.identifiability);
    end
end

% save('Seguin_Metrics_All.mat', 'results', '-v6');
% fprintf('\nAll Seguin metrics computed and saved to Seguin_Metrics_All.mat\n');

%% ===== Local helper (minimal variable-picker used above) =====
function M = pick_var_or_fallback(tmpStruct, candidate_varnames)
% Given a loaded struct tmpStruct and a cell of candidate varnames (strings),
% return the first matching variable. If none match, return the first numeric
% variable with ndims>=2 as a fallback. Throws error if nothing found.
M = [];
for k = 1:numel(candidate_varnames)
    vn = candidate_varnames{k};
    if isfield(tmpStruct, vn)
        M = tmpStruct.(vn);
        return;
    end
end
% fallback: first numeric ndims>=2
fns = fieldnames(tmpStruct);
for k = 1:numel(fns)
    cand = tmpStruct.(fns{k});
    if isnumeric(cand) && ndims(cand) >= 2
        M = cand;
        return;
    end
end
error('pick_var_or_fallback: no suitable numeric variable found in .mat file.');
end