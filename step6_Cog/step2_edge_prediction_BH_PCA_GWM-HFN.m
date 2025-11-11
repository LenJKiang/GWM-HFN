%% ========================================================================
%  Brain Basis Set (BBS) prediction using gretna_prediction
%  Networks: GWM-HFN (GM-WM-GM) and GM-GM
%  Behaviors: ONLY those showing PLS significance
%  Author: Zeqiang Linli
%  Date:   2025-09-03
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------------
% (1) SETTINGS
% -------------------------------------------------------------------------
base_path        = 'E:\Neuroimage\MyProject\GMWM_Network\Data\BGSP\';

% List of behavior variables (must match filenames in ./Behav_Data/)
behavior_list    = {'EstIQ_Shipley_Int_Bin','Shipley_Vocab_Raw'};
behavior_files   = strcat(behavior_list, '.csv');
beh_n            = numel(behavior_list);

% ROI / model settings
roi_num          = 90;
K_fold           = 10;
C_type           = 'PCA';
C_thr            = 80;          % PCA cumulative variance threshold (%)
R_type           = 'glm';
Ncomp            = [];
Lambda           = [];

% Repetition and permutation parameters
rep_time         = 100;         % Repeated CV times
rand_time        = 10000;       % Permutation count

% Edge importance thresholds
threshold_k      = 3;           % Z threshold for edge importance
perc_pos         = 99;          % Positive tail percentile
perc_neg         = 1;           % Negative tail percentile

% Reproducibility
rng_master       = 20250903;
rng(rng_master,'twister');

%% ------------------------------------------------------------------------
% (2) LOAD CONNECTOME DATA & BUILD EDGE MATRICES
% -------------------------------------------------------------------------
load('E:\Neuroimage\MyProject\GMWM_Network\Data\BGSP\GM_WM_GM_Connections_array.mat', 'GM_WM_GM_conn_array');   % [90 x 90 x N]
load('E:\Neuroimage\MyProject\GMWM_Network\Data\BGSP\GM_GM_Connections_array.mat',   'GM_GM_conn_array');       % [90 x 90 x N]
IDs_struct = load('E:\Neuroimage\MyProject\GMWM_Network\Data\BGSP\subjectIDs.mat','modified_subjects');
ID_all     = IDs_struct.modified_subjects;

if size(GM_WM_GM_conn_array,3) ~= size(GM_GM_conn_array,3)
    error('Subject count mismatch between GM-WM-GM and GM-GM arrays.');
end
nSub = size(GM_WM_GM_conn_array,3);

% Lower triangle mask (excluding diagonal)
Mask_net   = logical(tril(ones(roi_num), -1));
edge_count = sum(Mask_net(:));
[row_idx, col_idx] = find(Mask_net);

% Preallocate predictor matrices (subjects x edges)
Predictors_GWMHFN = zeros(nSub, edge_count);
Predictors_GMGM   = zeros(nSub, edge_count);

% Flatten lower triangle for each subject
for s = 1:nSub
    A1 = GM_WM_GM_conn_array(:,:,s);
    A2 = GM_GM_conn_array(:,:,s);
    Predictors_GWMHFN(s,:) = A1(Mask_net);
    Predictors_GMGM(s,:)   = A2(Mask_net);
end

% Standardize predictors across subjects (column-wise)
Predictors_GWMHFN = zscore(Predictors_GWMHFN);
Predictors_GMGM   = zscore(Predictors_GMGM);

%% ------------------------------------------------------------------------
% (3) PREALLOCATE RESULTS CONTAINER
% -------------------------------------------------------------------------
template_struct = struct( ...
    'Behavior','', ...
    'Network','', ...                      % 'GWM-HFN' or 'GM-GM'
    'nSubjects',[], ...
    'R_all',zeros(rep_time,1), ...
    'P_all',zeros(rep_time,1), ...
    'ExplainedVar_all',zeros(rep_time,1), ...
    'R_mean',[], 'ExplainedVar_mean',[], 'P_cv_mean',[], ...
    'R_perm_dist',zeros(rand_time,1), ...
    'P_perm',[], ...
    'ConsensusWeightsAll',[], ...          % (edges x rep_time)
    'ConsensusWeightsMean',[], ...         % (edges x 1)
    'Important_Z_Idx',[], 'Important_Z_Weights',[], ...
    'Important_Perc_Idx',[], 'Important_Perc_Weights',[], ...
    'Z_thr',[], 'Perc_Pos',[], 'Perc_Neg',[]);

results_BBS = repmat(template_struct, beh_n*2, 1);  % Two rows per behavior (one per network)

%% ------------------------------------------------------------------------
% (4) MAIN LOOP OVER BEHAVIORS AND NETWORKS
% -------------------------------------------------------------------------
fprintf('=== Starting BBS (gretna_prediction) for %d behaviors × 2 networks ===\n', beh_n);
tic_total = tic;

for b = 1:beh_n
    beh_name  = behavior_list{b};
    beh_file  = fullfile(base_path,'Behav_Data',behavior_files{b});
    T = readtable(beh_file,'Delimiter',',');

    % Subject alignment between connectivity and behavior
    [isMatch, idxMatch] = ismember(ID_all, T.Match_ID);
    Y_full = T{idxMatch(isMatch),2};
    if any(isnan(Y_full)), Y_full(isnan(Y_full)) = nanmean(Y_full); end
    n_use = numel(Y_full);

    % Extract matched predictor rows
    X_GWMHFN = Predictors_GWMHFN(isMatch,:);
    X_GMGM   = Predictors_GMGM(isMatch,:);

    seed_behavior = 100000 + b*1000;  % Behavior-specific base seed
    local_blocks = cell(1,2);         % To hold network results (GWM-HFN / GM-GM)

    for net_i = 1:2
        if net_i == 1
            net_name = 'GWM-HFN';
            X = X_GWMHFN;
        else
            net_name = 'GM-GM';
            X = X_GMGM;
        end

        % Preallocate arrays for repeated CV & permutation
        R_all         = zeros(rep_time,1);
        P_all         = zeros(rep_time,1);
        EV_all        = zeros(rep_time,1);
        Consensus_all = zeros(edge_count, rep_time);
        R_perm        = zeros(rand_time,1);

        % --- Repeated cross-validation predictions ---
        parfor rcv = 1:rep_time
            rng(seed_behavior + net_i*10 + rcv,'twister');
            Res = gretna_prediction(X, Y_full, K_fold, C_type, C_thr, R_type, Ncomp, Lambda);
            R_all(rcv)          = Res.R;
            P_all(rcv)          = Res.P;
            EV_all(rcv)         = mean(Res.Explained_variance);
            Consensus_all(:,rcv)= Res.Consensus_weights;
            if mod(rcv,ceil(rep_time/5))==0
                fprintf('[%s | %s] Repetition %d/%d completed.\n', beh_name, net_name, rcv, rep_time);
            end
        end

        % --- Permutation test (shuffle behavior vector) ---
        parfor pr = 1:rand_time
            rng(seed_behavior + net_i*100 + pr,'twister');
            idx_perm = randperm(n_use);
            if all(idx_perm==(1:n_use))
                idx_perm = randperm(n_use); % Extremely rare identity fallback
            end
            Y_perm = Y_full(idx_perm);
            ResP = gretna_prediction(X, Y_perm, K_fold, C_type, C_thr, R_type, Ncomp, Lambda);
            R_perm(pr) = ResP.R;
            if pr==1 || mod(pr,2000)==0
                fprintf('[%s | %s] Permutation %d/%d\n', beh_name, net_name, pr, rand_time);
            end
        end

        % --- Aggregate performance metrics ---
        R_mean  = mean(R_all);
        EV_mean = mean(EV_all);
        P_cv    = mean(P_all);
        P_perm  = (sum(R_perm >= R_mean)+1)/(rand_time+1);

        % --- Edge importance (mean consensus weights) ---
        mean_weights = mean(Consensus_all,2);
        z_weights    = zscore(mean_weights);

        % Z-threshold selection
        flag_z       = abs(z_weights) > threshold_k;
        imp_idx_z    = find(flag_z);
        imp_w_z      = mean_weights(imp_idx_z);

        % Percentile-based selection (separate positive / negative)
        pos_vals     = mean_weights(mean_weights>0);
        neg_vals     = mean_weights(mean_weights<0);
        if isempty(pos_vals), thr_pos = Inf; else, thr_pos = prctile(pos_vals, perc_pos); end
        if isempty(neg_vals), thr_neg = -Inf; else, thr_neg = prctile(neg_vals, perc_neg); end
        flag_perc    = (mean_weights >= thr_pos) | (mean_weights <= thr_neg);
        imp_idx_perc = find(flag_perc);
        imp_w_perc   = mean_weights(imp_idx_perc);

        % --- Populate local block ---
        blk = template_struct;
        blk.Behavior               = beh_name;
        blk.Network                = net_name;
        blk.nSubjects              = n_use;
        blk.R_all                  = R_all;
        blk.P_all                  = P_all;
        blk.ExplainedVar_all       = EV_all;
        blk.R_mean                 = R_mean;
        blk.ExplainedVar_mean      = EV_mean;
        blk.P_cv_mean              = P_cv;
        blk.R_perm_dist            = R_perm;
        blk.P_perm                 = P_perm;
        blk.ConsensusWeightsAll    = Consensus_all;
        blk.ConsensusWeightsMean   = mean_weights;
        blk.Important_Z_Idx        = imp_idx_z;
        blk.Important_Z_Weights    = imp_w_z;
        blk.Important_Perc_Idx     = imp_idx_perc;
        blk.Important_Perc_Weights = imp_w_perc;
        blk.Z_thr                  = threshold_k;
        blk.Perc_Pos               = perc_pos;
        blk.Perc_Neg               = perc_neg;

        local_blocks{net_i} = blk;
    end

    % Assign back into master results (two consecutive rows per behavior)
    results_BBS( (b-1)*2 + 1 ) = local_blocks{1};
    results_BBS( (b-1)*2 + 2 ) = local_blocks{2};
end

elapsed_min = toc(tic_total)/60;
fprintf('=== All behaviors completed in %.1f minutes ===\n', elapsed_min);

%% ------------------------------------------------------------------------
% (5) EXPORT EDGE IMPORTANCE TABLES
% -------------------------------------------------------------------------
out_dir = fullfile(base_path,'BBS_Output');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

summary_cells = {
    'Behavior','Network','nSubjects','R_mean','ExplainedVar_mean','P_cv_mean','P_perm', ...
    'nImportant_Z','nImportant_Perc'};

for k = 1:numel(results_BBS)
    Rk = results_BBS(k);

    % Summary row
    summary_cells(end+1,:) = {
        Rk.Behavior, Rk.Network, Rk.nSubjects, Rk.R_mean, Rk.ExplainedVar_mean, ...
        Rk.P_cv_mean, Rk.P_perm, numel(Rk.Important_Z_Idx), numel(Rk.Important_Perc_Idx)
    };

    % Z-threshold edges
    if ~isempty(Rk.Important_Z_Idx)
        nodes1_z = row_idx(Rk.Important_Z_Idx);
        nodes2_z = col_idx(Rk.Important_Z_Idx);
        Tz = table(nodes1_z, nodes2_z, Rk.Important_Z_Weights, ...
            'VariableNames',{'Node1','Node2','Weight'});
        fn_z = fullfile(out_dir, sprintf('EdgeImportance_%s_%s_ZThr.csv', ...
            Rk.Behavior, strrep(Rk.Network,'-','')));
        writetable(Tz, fn_z);
    end

    % Percentile edges
    if ~isempty(Rk.Important_Perc_Idx)
        nodes1_p = row_idx(Rk.Important_Perc_Idx);
        nodes2_p = col_idx(Rk.Important_Perc_Idx);
        Tp = table(nodes1_p, nodes2_p, Rk.Important_Perc_Weights, ...
            'VariableNames',{'Node1','Node2','Weight'});
        fn_p = fullfile(out_dir, sprintf('EdgeImportance_%s_%s_Perc.csv', ...
            Rk.Behavior, strrep(Rk.Network,'-','')));
        writetable(Tp, fn_p);
    end
end

SummaryTable = cell2table(summary_cells(2:end,:), 'VariableNames', summary_cells(1,:));
writetable(SummaryTable, fullfile(out_dir,'BBS_Summary.csv'));

%% ------------------------------------------------------------------------
% (6) SAVE WORKSPACE RESULTS
% -------------------------------------------------------------------------
save(fullfile(out_dir,'BBS_Results_GWMHFN_GMGM.mat'), ...
    'results_BBS','behavior_list','Mask_net','row_idx','col_idx', ...
    'threshold_k','perc_pos','perc_neg','rep_time','rand_time','C_thr','K_fold','-v7.3');

fprintf('\nOutputs saved to: %s\n', out_dir);
fprintf('  - BBS_Results_GWMHFN_GMGM.mat\n  - BBS_Summary.csv\n  - EdgeImportance_* CSV files\n');

%% ------------------------------------------------------------------------
% (6b) EXPLICIT EXPORTS + INTEGRATED PERMUTATION FIGURE & INTERSECTION
%      (Translated & Modularization Enhancements Applied Above)
%      Retaining user-defined logic without behavioral change.
% -------------------------------------------------------------------------
% Parameters
behavior_vocab = 'Shipley_Vocab_Raw';
behavior_estIQ = 'EstIQ_Shipley_Int_Bin';
z_thr_export   = threshold_k;
out_bbs_simple = fullfile(out_dir,'BBS_SimpleExports');
if ~exist(out_bbs_simple,'dir'), mkdir(out_bbs_simple); end

% Locate entries
idx_vocab_gwmhfn = find(strcmp({results_BBS.Behavior},behavior_vocab) & strcmp({results_BBS.Network},'GWM-HFN'),1);
idx_vocab_gmgm   = find(strcmp({results_BBS.Behavior},behavior_vocab) & strcmp({results_BBS.Network},'GM-GM'),1);

% Mean consensus weights
weights_vocab_gwmhfn = results_BBS(idx_vocab_gwmhfn).ConsensusWeightsMean;
z_vocab_gwmhfn       = zscore(weights_vocab_gwmhfn);
flag_vocab_gwmhfn    = abs(z_vocab_gwmhfn) > z_thr_export;
sel_idx_vocab_gwmhfn = find(flag_vocab_gwmhfn);
sel_w_vocab_gwmhfn   = weights_vocab_gwmhfn(sel_idx_vocab_gwmhfn);

weights_vocab_gmgm   = results_BBS(idx_vocab_gmgm).ConsensusWeightsMean;
z_vocab_gmgm         = zscore(weights_vocab_gmgm);
flag_vocab_gmgm      = abs(z_vocab_gmgm) > z_thr_export;
sel_idx_vocab_gmgm   = find(flag_vocab_gmgm);
sel_w_vocab_gmgm     = weights_vocab_gmgm(sel_idx_vocab_gmgm);

% Build BrainNet edge matrices
E_vocab_gwmhfn = zeros(roi_num); 
r_sel_vocab_gwmhfn = row_idx(sel_idx_vocab_gwmhfn);
c_sel_vocab_gwmhfn = col_idx(sel_idx_vocab_gwmhfn);
lin_gwmhfn = sub2ind([roi_num roi_num], r_sel_vocab_gwmhfn, c_sel_vocab_gwmhfn);
E_vocab_gwmhfn(lin_gwmhfn) = sel_w_vocab_gwmhfn;
E_vocab_gwmhfn(sub2ind([roi_num roi_num], c_sel_vocab_gwmhfn, r_sel_vocab_gwmhfn)) = sel_w_vocab_gwmhfn;

E_vocab_gmgm = zeros(roi_num);
r_sel_vocab_gmgm = row_idx(sel_idx_vocab_gmgm);
c_sel_vocab_gmgm = col_idx(sel_idx_vocab_gmgm);
lin_gmgm = sub2ind([roi_num roi_num], r_sel_vocab_gmgm, c_sel_vocab_gmgm);
E_vocab_gmgm(lin_gmgm) = sel_w_vocab_gmgm;
E_vocab_gmgm(sub2ind([roi_num roi_num], c_sel_vocab_gmgm, r_sel_vocab_gmgm)) = sel_w_vocab_gmgm;

% Load AAL node file
node_table_simple = readtable('AAL90.node','FileType','text','Delimiter','\t','ReadVariableNames',false);
node_coord_simple = table2array(node_table_simple(:,1:5));
node_label_simple = table2cell(node_table_simple(:,6));

% Node degree -> node size mapping
deg_gwmhfn = sum(E_vocab_gwmhfn~=0,2);
deg_gmgm   = sum(E_vocab_gmgm~=0,2);
max_deg_gwmhfn = max(deg_gwmhfn);
max_deg_gmgm   = max(deg_gmgm);
node_size_gwmhfn = zeros(roi_num,1);
node_size_gmgm   = zeros(roi_num,1);
if max_deg_gwmhfn>0
    node_size_gwmhfn = 5 * deg_gwmhfn / max_deg_gwmhfn;
end
if max_deg_gmgm>0
    node_size_gmgm   = 5 * deg_gmgm   / max_deg_gmgm;
end

% Write BrainNet files (GWM-HFN)
edge_file_gwmhfn = fullfile(out_bbs_simple, sprintf('%s_GWMHFN_Zgt%d.edge',behavior_vocab,z_thr_export));
dlmwrite(edge_file_gwmhfn,E_vocab_gwmhfn,'delimiter','\t','precision','%.6f');
node_file_gwmhfn = strrep(edge_file_gwmhfn,'.edge','.node');
fid_gwmhfn = fopen(node_file_gwmhfn,'w');
for nn = 1:roi_num
    fprintf(fid_gwmhfn,'%.4f\t%.4f\t%.4f\t%d\t%.4f\t%s\n', ...
        node_coord_simple(nn,1), node_coord_simple(nn,2), node_coord_simple(nn,3), ...
        node_coord_simple(nn,4), node_size_gwmhfn(nn), node_label_simple{nn});
end
fclose(fid_gwmhfn);

% Write BrainNet files (GM-GM)
edge_file_gmgm = fullfile(out_bbs_simple, sprintf('%s_GMGM_Zgt%d.edge',behavior_vocab,z_thr_export));
dlmwrite(edge_file_gmgm,E_vocab_gmgm,'delimiter','\t','precision','%.6f');
node_file_gmgm = strrep(edge_file_gmgm,'.edge','.node');
fid_gmgm = fopen(node_file_gmgm,'w');
for nn = 1:roi_num
    fprintf(fid_gmgm,'%.4f\t%.4f\t%.4f\t%d\t%.4f\t%s\n', ...
        node_coord_simple(nn,1), node_coord_simple(nn,2), node_coord_simple(nn,3), ...
        node_coord_simple(nn,4), node_size_gmgm(nn), node_label_simple{nn});
end
fclose(fid_gmgm);

% Permutation distributions (already computed) - combine figure
perm_gwmhfn = results_BBS(idx_vocab_gwmhfn).R_perm_dist;
obs_gwmhfn  = results_BBS(idx_vocab_gwmhfn).R_mean;
p_gwm       = results_BBS(idx_vocab_gwmhfn).P_perm;
perm_sorted_gwmhfn = sort(perm_gwmhfn);
crit95_gwmhfn = perm_sorted_gwmhfn(round(0.95*length(perm_sorted_gwmhfn)));

perm_gmgm = results_BBS(idx_vocab_gmgm).R_perm_dist;
obs_gmgm  = results_BBS(idx_vocab_gmgm).R_mean;
p_gmgm    = results_BBS(idx_vocab_gmgm).P_perm;
perm_sorted_gmgm = sort(perm_gmgm);
crit95_gmgm = perm_sorted_gmgm(round(0.95*length(perm_sorted_gmgm)));

% Density (auto bandwidth)
[f1,x1] = ksdensity(perm_gwmhfn);
[f2,x2] = ksdensity(perm_gmgm);
sig_mask1 = x1 >= crit95_gwmhfn;
sig_mask2 = x2 >= crit95_gmgm;
ymax = max([f1(:); f2(:)]) * 1.20;

fig_comb = figure('Color','white');
set(fig_comb,'Units','pixels','Position',[200 200 800 420]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
area(x1,f1,'FaceColor',[0.78 0.88 1],'EdgeColor','none','FaceAlpha',0.65); hold on;
area(x1(sig_mask1),f1(sig_mask1),'FaceColor',[0.45 0.70 1],'EdgeColor','none','FaceAlpha',0.70);
plot([obs_gwmhfn obs_gwmhfn],[0 ymax*0.92],'r-','LineWidth',2);
set(gca,'LineWidth',1.2,'FontSize',12,'TickDir','out','Box','on');
xlabel('Correlation Coefficient (R)','FontSize',15,'FontWeight','bold');
ylabel('Density','FontSize',15,'FontWeight','bold');
title('Permutation (GWM-HFN)','FontSize',17,'FontWeight','bold');
ylim([0 ymax]);
legend({'Distribution','>95% Region',sprintf('Observed R = %.3f\n(p = %.4f)',obs_gwmhfn,p_gwm)}, ...
       'Location','NorthWest','Box','off','FontSize',12,'FontWeight','bold');

nexttile;
area(x2,f2,'FaceColor',[0.90 0.85 1],'EdgeColor','none','FaceAlpha',0.65); hold on;
area(x2(sig_mask2),f2(sig_mask2),'FaceColor',[0.65 0.55 0.95],'EdgeColor','none','FaceAlpha',0.70);
plot([obs_gmgm obs_gmgm],[0 ymax*0.92],'r-','LineWidth',2);
set(gca,'LineWidth',1.2,'FontSize',12,'TickDir','out','Box','on');
xlabel('Correlation Coefficient (R)','FontSize',15,'FontWeight','bold');
ylabel('Density','FontSize',15,'FontWeight','bold');
title('Permutation (GM-GM)','FontSize',17,'FontWeight','bold');
ylim([0 ymax]);
legend({'Distribution','>95% Region',sprintf('Observed R = %.3f\n(p = %.4f)',obs_gmgm,p_gmgm)}, ...
       'Location','NorthWest','Box','off','FontSize',12,'FontWeight','bold');

out_comb_png = fullfile(out_bbs_simple, sprintf('PermDist_%s_GWMHFN_GMGM_Zgt%d.png', behavior_vocab, z_thr_export));
out_comb_tif = strrep(out_comb_png,'.png','.tiff');
print(fig_comb, out_comb_png,'-dpng','-r300');
print(fig_comb, out_comb_tif,'-dtiff','-r300');

PermSource = table( ...
    [x1(:); NaN; x2(:)], ...
    [f1(:); NaN; f2(:)], ...
    [repmat({'GWM-HFN'},numel(x1),1); {'---'}; repmat({'GM-GM'},numel(x2),1)], ...
    'VariableNames',{'X','Density','Network'});
writetable(PermSource, fullfile(out_bbs_simple, sprintf('PermDist_%s_GWMHFN_GMGM_Zgt%d_SourceData.csv', behavior_vocab, z_thr_export)));

StatsSummary = table( ...
    {'GWM-HFN';'GM-GM'}, ...
    [obs_gwmhfn; obs_gmgm], ...
    [p_gwm; p_gmgm], ...
    [crit95_gwmhfn; crit95_gmgm], ...
    'VariableNames',{'Network','Observed_R','Perm_P','CritR_95'});
writetable(StatsSummary, fullfile(out_bbs_simple, sprintf('PermDist_%s_GWMHFN_GMGM_Zgt%d_Stats.csv', behavior_vocab, z_thr_export)));

% Framework-level intersection for the same behavior
edges_vocab_gwmhfn = sel_idx_vocab_gwmhfn;
edges_vocab_gmgm   = sel_idx_vocab_gmgm;
edges_inter_framework = intersect(edges_vocab_gwmhfn, edges_vocab_gmgm);

E_inter_framework = zeros(roi_num);
r_int_fw = row_idx(edges_inter_framework);
c_int_fw = col_idx(edges_inter_framework);
lin_int_fw = sub2ind([roi_num roi_num], r_int_fw, c_int_fw);
w_gwmhfn_int = weights_vocab_gwmhfn(edges_inter_framework);
w_gmgm_int   = weights_vocab_gmgm(edges_inter_framework);
w_mean_int_fw = 0.5 * (w_gwmhfn_int + w_gmgm_int);

E_inter_framework(lin_int_fw) = w_mean_int_fw;
E_inter_framework(sub2ind([roi_num roi_num], c_int_fw, r_int_fw)) = w_mean_int_fw;

edge_file_inter_fw = fullfile(out_bbs_simple, ...
    sprintf('Intersection_%s_GWMHFN_GMGM_Zgt%d.edge', behavior_vocab, z_thr_export));
dlmwrite(edge_file_inter_fw, E_inter_framework, 'delimiter','\t','precision','%.6f');

deg_int_fw = sum(E_inter_framework~=0,2);
max_deg_int_fw = max(deg_int_fw);
node_size_int_fw = zeros(roi_num,1);
if max_deg_int_fw>0
    node_size_int_fw = 5 * deg_int_fw / max_deg_int_fw;
end
node_file_inter_fw = strrep(edge_file_inter_fw,'.edge','.node');
fid_ifw = fopen(node_file_inter_fw,'w');
for nn = 1:roi_num
    fprintf(fid_ifw,'%.4f\t%.4f\t%.4f\t%d\t%.4f\t%s\n', ...
        node_coord_simple(nn,1), node_coord_simple(nn,2), node_coord_simple(nn,3), ...
        node_coord_simple(nn,4), node_size_int_fw(nn), node_label_simple{nn});
end
fclose(fid_ifw);

T_inter_fw = table( ...
    r_int_fw(:), c_int_fw(:), ...
    w_gwmhfn_int(:), w_gmgm_int(:), w_mean_int_fw(:), ...
    'VariableNames', {'ROI_i','ROI_j','Weight_GWMHFN','Weight_GMGM','MeanWeight'});
csv_inter_fw = fullfile(out_bbs_simple, ...
    sprintf('IntersectionEdges_%s_GWMHFN_GMGM_Zgt%d.csv', behavior_vocab, z_thr_export));
writetable(T_inter_fw, csv_inter_fw);

% Circos link for intersection
genes_file_circos = 'E:\Neuroimage\Software\circos-0.69-6\AAL\genes.labelsAAL.links.txt';
genesInfo_circos = readtable(genes_file_circos,'FileType','text','Delimiter','\t','ReadVariableNames',false);
genesInfo_circos.Properties.VariableNames = {'Network','Start','End','Region'};

w_vec_fw = w_mean_int_fw(:);
r_vec_fw = r_int_fw(:);
c_vec_fw = c_int_fw(:);
den_fw = max(abs(w_vec_fw)); den_fw = den_fw + (den_fw==0);
thick_fw = 1 + 9 * abs(w_vec_fw) / den_fw;
pos_mask_fw = (w_vec_fw > 0);
colorR_fw = 255 * pos_mask_fw;
colorG_fw = zeros(size(w_vec_fw));
colorB_fw = 255 * (~pos_mask_fw);
link_file_inter_fw = fullfile(out_bbs_simple, ...
    sprintf('links_Intersection_%s_GWMHFN_GMGM_Zgt%d.txt', behavior_vocab, z_thr_export));
fid_lfw = fopen(link_file_inter_fw,'w');
for ii = 1:numel(w_vec_fw)
    gA = genesInfo_circos(r_vec_fw(ii),:);
    gB = genesInfo_circos(c_vec_fw(ii),:);
    fprintf(fid_lfw,'%s %d %d %s %d %d color=%d,%d,%d,thickness=%.2f\n', ...
        gA.Network{1}, gA.Start, gA.End, gB.Network{1}, gB.Start, gB.End, ...
        colorR_fw(ii), colorG_fw(ii), colorB_fw(ii), thick_fw(ii));
end
fclose(fid_lfw);

J_framework = numel(edges_inter_framework) / numel(union(edges_vocab_gwmhfn, edges_vocab_gmgm));
D_framework = 2 * numel(edges_inter_framework) / (numel(edges_vocab_gwmhfn) + numel(edges_vocab_gmgm));
fprintf('  Jaccard=%.4f  Dice=%.4f\n', J_framework, D_framework);

% Individual Circos link files for single frameworks (kept as in prior logic)
% GWM-HFN
w_vec_1 = sel_w_vocab_gwmhfn(:);
r_vec_1 = r_sel_vocab_gwmhfn(:);
c_vec_1 = c_sel_vocab_gwmhfn(:);
den_1   = max(abs(w_vec_1)); den_1 = den_1 + (den_1==0);
thick_1 = 1 + 9 * abs(w_vec_1) / den_1;
pos_mask_1 = (w_vec_1 > 0);
colorR_1 = 255 * pos_mask_1;
colorG_1 = zeros(size(w_vec_1));
colorB_1 = 255 * (~pos_mask_1);
link_file_1 = fullfile(out_bbs_simple, sprintf('links_%s_GWMHFN_Zgt%d.txt', behavior_vocab, z_thr_export));
fid_1 = fopen(link_file_1,'w');
for ii = 1:numel(w_vec_1)
    gA = genesInfo_circos(r_vec_1(ii),:);
    gB = genesInfo_circos(c_vec_1(ii),:);
    fprintf(fid_1,'%s %d %d %s %d %d color=%d,%d,%d,thickness=%.2f\n', ...
        gA.Network{1}, gA.Start, gA.End, gB.Network{1}, gB.Start, gB.End, ...
        colorR_1(ii), colorG_1(ii), colorB_1(ii), thick_1(ii));
end
fclose(fid_1);

% GM-GM
w_vec_2 = sel_w_vocab_gmgm(:);
r_vec_2 = r_sel_vocab_gmgm(:);
c_vec_2 = c_sel_vocab_gmgm(:);
den_2   = max(abs(w_vec_2)); den_2 = den_2 + (den_2==0);
thick_2 = 1 + 9 * abs(w_vec_2) / den_2;
pos_mask_2 = (w_vec_2 > 0);
colorR_2 = 255 * pos_mask_2;
colorG_2 = zeros(size(w_vec_2));
colorB_2 = 255 * (~pos_mask_2);
link_file_2 = fullfile(out_bbs_simple, sprintf('links_%s_GMGM_Zgt%d.txt', behavior_vocab, z_thr_export));
fid_2 = fopen(link_file_2,'w');
for ii = 1:numel(w_vec_2)
    gA = genesInfo_circos(r_vec_2(ii),:);
    gB = genesInfo_circos(c_vec_2(ii),:);
    fprintf(fid_2,'%s %d %d %s %d %d color=%d,%d,%d,thickness=%.2f\n', ...
        gA.Network{1}, gA.Start, gA.End, gB.Network{1}, gB.Start, gB.End, ...
        colorR_2(ii), colorG_2(ii), colorB_2(ii), thick_2(ii));
end
fclose(fid_2);

fprintf('[Combined Permutation Figure] saved:\n  %s\n  %s\n', out_comb_png, out_comb_tif);
fprintf('[Intersection Framework] %s: %d common high-weight edges (Z>|%g|)\n', ...
    behavior_vocab, numel(edges_inter_framework), z_thr_export);
fprintf('  BrainNet edge: %s\n', edge_file_inter_fw);
fprintf('  BrainNet node: %s\n', node_file_inter_fw);
fprintf('  Edge list CSV: %s\n', csv_inter_fw);
fprintf('  Circos link (intersection): %s\n', link_file_inter_fw);

fprintf('[Circos Links] %s\n', link_file_1);
fprintf('[Circos Links] %s\n', link_file_2);

fprintf('\n[Simple Export] %s (GWM-HFN): %d edges (Z>|%g|)\n',behavior_vocab,numel(sel_idx_vocab_gwmhfn),z_thr_export);
fprintf('[Simple Export] %s (GM-GM):   %d edges (Z>|%g|)\n',behavior_vocab,numel(sel_idx_vocab_gmgm),z_thr_export);

%% ------------------------------------------------------------------------
% (7) QUICK TEXT SUMMARY
% -------------------------------------------------------------------------
fprintf('\n=== Quick Performance Overview ===\n');
for k = 1:numel(results_BBS)
    Rk = results_BBS(k);
    fprintf('%-25s | %-7s | R=%.3f (perm p=%.4f) EV=%.3f ZEdges=%d PercEdges=%d\n', ...
        Rk.Behavior, Rk.Network, Rk.R_mean, Rk.P_perm, ...
        Rk.ExplainedVar_mean, numel(Rk.Important_Z_Idx), numel(Rk.Important_Perc_Idx));
end
fprintf('=================================================================\n');

%% ------------------------------------------------------------------------
