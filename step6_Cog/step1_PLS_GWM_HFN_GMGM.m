%% ========================================================================
%  Step 6. Cognitive Association Analysis 
%  Author: Zeqiang Linli  |  Date: 2025-08-21
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------------
%  (0) Paths and basic settings
% -------------------------------------------------------------------------
cd 'E:\Neuroimage\MyProject\GMWM_Network\Code_R1\step6_Cog'
datapath        = 'E:\Neuroimage\MyProject\GMWM_Network\Data\BGSP\';
roi_num         = 90;
rand_time       = 10000;
do_bootstrap    = true;
B_boot          = 1000;
store_full_boot = false;
base_seed_perm  = 120000;
base_seed_boot  = 220000;

% Behavioral files
BH_name  = {dir(fullfile(datapath,'Behav_Data','*.csv')).name};
BH_index = 1:numel(BH_name);
fprintf('Found %d behavioral variables.\n', numel(BH_name));

% Lower triangle mask
Mask_net = logical(tril(ones(roi_num), -1));
nEdges   = sum(Mask_net(:));
[r_all, c_all] = find(Mask_net);

% Parallel pool
if isempty(gcp('nocreate')), parpool('local'); end

%% ------------------------------------------------------------------------
%  Load BrainNet node template
% -------------------------------------------------------------------------
node_template_file = 'AAL90.node';
tbl_node = readtable(node_template_file,'FileType','text','Delimiter','\t','ReadVariableNames',false);
node_base   = table2array(tbl_node(:,1:5));
node_labels = table2cell(tbl_node(:,6));
if size(node_base,1)~=roi_num
    error('Node template count mismatch.');
end
MaxNodeSize = 5;

%% ========================================================================
%  PART A. GWM-HFN PLS
% ========================================================================
fprintf('\n================= PART A: GWM-HFN (GM-WM-GM) PLS =================\n');
load(fullfile(datapath,'GM_WM_GM_Connections_array.mat'));
ID_struct = load(fullfile(datapath,'subjectIDs.mat'));
ID = ID_struct.modified_subjects;
nsubjs = size(GM_WM_GM_conn_array,3);
fprintf('Loaded GM-WM-GM connectivity: %d subjects.\n', nsubjs);

Mat_data = zeros(nsubjs, nEdges);
for s = 1:nsubjs
    M = GM_WM_GM_conn_array(:,:,s);
    Mat_data(s,:) = M(Mask_net);
end
Mat_data = zscore(Mat_data);

results = repmat(struct( ...
    'BH_name','', 'nsubj',[], 'XS',[], 'YS',[], ...
    'variance_X',[], 'variance_Y',[], 'weights',[], ...
    'variance_Y_rand',zeros(rand_time,1), 'p',[], ...
    'boot_n',[], 'boot_mean',[], 'boot_SE',[], 'boot_BR',[], ...
    'boot_signProp',[], 'boot_flag_BR3',[], 'boot_p_norm',[], 'boot_p_FDR',[]), numel(BH_index),1);

parfor i = BH_index
    rng(base_seed_perm + i,'twister');
    bh_file = fullfile(datapath,'Behav_Data',BH_name{i});
    BH_Data = readtable(bh_file,'Delimiter',',');

    [isMatch, idxMatch] = ismember(ID, BH_Data.Match_ID);
    X = Mat_data(isMatch,:);
    Y = BH_Data{idxMatch(isMatch),2};
    X = zscore(X);
    Y(isnan(Y)) = nanmean(Y);
    nS = length(Y);

    [~,~,XS,YS,~,PCTVAR,~,stats] = plsregress(X,zscore(Y),1);

    perm_varY = zeros(rand_time,1);
    for ip = 1:rand_time
        rp = randperm(nS);
        if all(rp==(1:nS)), rp = randperm(nS); end
        Yp = Y(rp);
        [~,~,~,~,~,PCTVAR_p,~,~] = plsregress(X,zscore(Yp),1);
        perm_varY(ip) = PCTVAR_p(2,:);
    end

    boot_struct = struct('boot_n',[],'boot_mean',[],'boot_SE',[],'boot_BR',[], ...
        'boot_signProp',[],'boot_flag_BR3',[],'boot_p_norm',[],'boot_p_FDR',[]);
    if do_bootstrap
        rng(base_seed_boot + i,'twister');
        boot_out = bootstrap_pls_weights(X,Y,stats.W,B_boot,store_full_boot);
        fn = fieldnames(boot_out);
        for f = 1:numel(fn)
            boot_struct.(fn{f}) = boot_out.(fn{f});
        end
    end

    rtmp = results(i);
    rtmp.BH_name          = erase(BH_name{i},'.csv');
    rtmp.nsubj            = nS;
    rtmp.XS               = XS;
    rtmp.YS               = YS;
    rtmp.variance_X       = PCTVAR(1,:);
    rtmp.variance_Y       = PCTVAR(2,:);
    rtmp.weights          = stats.W;
    rtmp.variance_Y_rand  = perm_varY;
    rtmp.p                = (sum(perm_varY >= rtmp.variance_Y)+1)/(rand_time+1);
    rtmp.boot_n           = boot_struct.boot_n;
    rtmp.boot_mean        = boot_struct.boot_mean;
    rtmp.boot_SE          = boot_struct.boot_SE;
    rtmp.boot_BR          = boot_struct.boot_BR;
    rtmp.boot_signProp    = boot_struct.boot_signProp;
    rtmp.boot_flag_BR3    = boot_struct.boot_flag_BR3;
    rtmp.boot_p_norm      = boot_struct.boot_p_norm;
    rtmp.boot_p_FDR       = boot_struct.boot_p_FDR;

    results(i) = rtmp;
    fprintf('[GWM-HFN] %d/%d %s: VarY=%.4f p=%.4g Stable(|BR|>3)=%d\n', ...
        i, numel(BH_name), rtmp.BH_name, rtmp.variance_Y, rtmp.p, sum(rtmp.boot_flag_BR3==1));
end

stable_table_A = table();
for i = BH_index
    if ~isempty(results(i).boot_BR)
        idx_stable = find(results(i).boot_flag_BR3);
        if ~isempty(idx_stable)
            addT = table(repmat({results(i).BH_name},numel(idx_stable),1), idx_stable, ...
                results(i).weights(idx_stable), results(i).boot_SE(idx_stable), ...
                results(i).boot_BR(idx_stable), results(i).boot_signProp(idx_stable), ...
                results(i).boot_p_FDR(idx_stable), ...
                'VariableNames',{'Behavior','EdgeIndex','Weight','SE','BR','SignStability','FDRp'});
            stable_table_A = [stable_table_A; addT]; %#ok<AGROW>
        end
    end
end
writetable(stable_table_A,'PLS_Bootstrap_StableEdges_GWMHFN.csv');
p_fdr_GWMHFN = mafdr([results.p],'BHFDR',true);
save('Edge_BH_PLS_GMWMGM_BGSP.mat','results','p_fdr_GWMHFN');

%% ========================================================================
%  PART B. GM-GM PLS
% ========================================================================
fprintf('\n================= PART B: GM-GM PLS (Benchmark) ===================\n');
load(fullfile(datapath,'GM_GM_Connections_array.mat'));
nsubjs_b = size(GM_GM_conn_array,3);
fprintf('Loaded GM-GM connectivity: %d subjects.\n', nsubjs_b);

Mat_data_B = zeros(nsubjs_b, nEdges);
for s = 1:nsubjs_b
    M = GM_GM_conn_array(:,:,s);
    Mat_data_B(s,:) = M(Mask_net);
end
Mat_data_B = zscore(Mat_data_B);

results_GMGM = repmat(struct( ...
    'BH_name','', 'nsubj',[], 'XS',[], 'YS',[], ...
    'variance_X',[], 'variance_Y',[], 'weights',[], ...
    'variance_Y_rand',zeros(rand_time,1), 'p',[], ...
    'boot_n',[], 'boot_mean',[], 'boot_SE',[], 'boot_BR',[], ...
    'boot_signProp',[], 'boot_flag_BR3',[], 'boot_p_norm',[], 'boot_p_FDR',[]), numel(BH_index),1);

parfor i = BH_index
    rng(base_seed_perm + 5000 + i,'twister');
    bh_file = fullfile(datapath,'Behav_Data',BH_name{i});
    BH_Data = readtable(bh_file,'Delimiter',',');

    [isMatch, idxMatch] = ismember(ID, BH_Data.Match_ID);
    X = Mat_data_B(isMatch,:);
    Y = BH_Data{idxMatch(isMatch),2};
    X = zscore(X);
    Y(isnan(Y)) = nanmean(Y);
    nS = length(Y);

    [~,~,XS,YS,~,PCTVAR,~,stats] = plsregress(X,zscore(Y),1);

    perm_varY = zeros(rand_time,1);
    for ip = 1:rand_time
        rp = randperm(nS);
        if all(rp==(1:nS)), rp = randperm(nS); end
        Yp = Y(rp);
        [~,~,~,~,~,PCTVAR_p,~,~] = plsregress(X,zscore(Yp),1);
        perm_varY(ip) = PCTVAR_p(2,:);
    end

    boot_struct = struct('boot_n',[],'boot_mean',[],'boot_SE',[],'boot_BR',[], ...
        'boot_signProp',[],'boot_flag_BR3',[],'boot_p_norm',[],'boot_p_FDR',[]);
    if do_bootstrap
        rng(base_seed_boot + 5000 + i,'twister');
        boot_out = bootstrap_pls_weights(X,Y,stats.W,B_boot,store_full_boot);
        fn = fieldnames(boot_out);
        for f = 1:numel(fn)
            boot_struct.(fn{f}) = boot_out.(fn{f});
        end
    end

    rtmp = results_GMGM(i);
    rtmp.BH_name          = erase(BH_name{i},'.csv');
    rtmp.nsubj            = nS;
    rtmp.XS               = XS;
    rtmp.YS               = YS;
    rtmp.variance_X       = PCTVAR(1,:);
    rtmp.variance_Y       = PCTVAR(2,:);
    rtmp.weights          = stats.W;
    rtmp.variance_Y_rand  = perm_varY;
    rtmp.p                = (sum(perm_varY >= rtmp.variance_Y)+1)/(rand_time+1);
    rtmp.boot_n           = boot_struct.boot_n;
    rtmp.boot_mean        = boot_struct.boot_mean;
    rtmp.boot_SE          = boot_struct.boot_SE;
    rtmp.boot_BR          = boot_struct.boot_BR;
    rtmp.boot_signProp    = boot_struct.boot_signProp;
    rtmp.boot_flag_BR3    = boot_struct.boot_flag_BR3;
    rtmp.boot_p_norm      = boot_struct.boot_p_norm;
    rtmp.boot_p_FDR       = boot_struct.boot_p_FDR;

    results_GMGM(i) = rtmp;
    fprintf('[GM-GM] %d/%d %s: VarY=%.4f p=%.4g Stable(|BR|>3)=%d\n', ...
        i, numel(BH_name), rtmp.BH_name, rtmp.variance_Y, rtmp.p, sum(rtmp.boot_flag_BR3==1));
end

stable_table_B = table();
for i = BH_index
    if ~isempty(results_GMGM(i).boot_BR)
        idx_stable = find(results_GMGM(i).boot_flag_BR3);
        if ~isempty(idx_stable)
            addT = table(repmat({results_GMGM(i).BH_name},numel(idx_stable),1), idx_stable, ...
                results_GMGM(i).weights(idx_stable), results_GMGM(i).boot_SE(idx_stable), ...
                results_GMGM(i).boot_BR(idx_stable), results_GMGM(i).boot_signProp(idx_stable), ...
                results_GMGM(i).boot_p_FDR(idx_stable), ...
                'VariableNames',{'Behavior','EdgeIndex','Weight','SE','BR','SignStability','FDRp'});
            stable_table_B = [stable_table_B; addT]; %#ok<AGROW>
        end
    end
end
writetable(stable_table_B,'PLS_Bootstrap_StableEdges_GMGM.csv');
p_fdr_GMGM = mafdr([results_GMGM.p],'BHFDR',true);
save('Edge_BH_PLS_GMGM_BGSP.mat','results_GMGM','p_fdr_GMGM');

%% ========================================================================
%  PART C. Scatterplot visualizations (retain original style)
% ========================================================================
load('Edge_BH_PLS_GMWMGM_BGSP.mat');  % results (GWM-HFN)
S = load('../../Data/BGSP/subjectIDs.mat'); 
ID = S.modified_subjects;

target_names = {'EstIQ_Shipley_Int_Bin','Shipley_Vocab_Raw','EstIQ_Matrix_Int_Bin','Matrix_WAIS'};
all_names = {results.BH_name};
sel_mask = ismember(all_names, target_names);
filtered_results_GWM = results(sel_mask);

% Prepare Excel for scatter plot source data
scatter_xlsx = 'PLS_Scatter_SourceData.xlsx';
if exist(scatter_xlsx,'file'), delete(scatter_xlsx); end

% Figure 1: Shipley_Vocab_Raw (GWM-HFN)
selected_row = filtered_results_GWM(4);
X = selected_row.XS;
T = readtable(fullfile('../../Data/BGSP/','Behav_Data',[target_names{2} '.csv']),'Delimiter',',');
[isMatch, matchIdx] = ismember(ID, T.Match_ID);
Y = T{matchIdx(isMatch),2};
figure; out = scatplot(X, Y, [], [], [], [], [], 22); %#ok<NASGU>
colorbar('off'); hold on;
mdl = fitlm(X, Y);
x_range = linspace(min(X), max(X), 100);
y_fit = predict(mdl, x_range');
plot(x_range, y_fit, 'LineWidth', 2.5, 'Color', 'red');
[r_ship_gwm, p_ship_gwm] = corr(Y, X, 'type','Pearson');
text(0.05, 0.95, sprintf('r = %.3f', r_ship_gwm), 'Units','normalized','FontSize',20,'Color','black','FontWeight','bold');
text(0.05, 0.90, sprintf('p = %.3g', p_ship_gwm), 'Units','normalized','FontSize',20,'Color','black','FontWeight','bold');
title('Correlation of PLS1 score and Shipley vocabulary task', ...
    'FontSize',23,'FontWeight','bold','Units','normalized','Position',[0.5, 1.02]);
xlabel('PLS1 score (GWM-HFN)','FontSize',21,'FontWeight','bold');
ylabel('EstIQ (Shipley)','FontSize',21,'FontWeight','bold');
xlim([-0.15, 0.12]); ylim([20, 45]);
xticks(-0.15:0.05:0.12); yticks(20:5:45);
set(gca,'FontSize',16,'FontWeight','bold'); grid on;
set(gcf,'Position',[1200, 600, 700, 600]);
print(fullfile('./Figures/','Corr_PLS1_Shipley_Vocab_Raw_GWM-HFN.png'),'-dpng','-r300'); hold off;

% === NEW EXPORT === Scatter source data (Sheet 1)
tbl_scatter1 = table(ID(isMatch), X, Y, 'VariableNames',{'SubjectID','PLS1_Score','Shipley_Vocab_Raw'});
tbl_meta1 = table(r_ship_gwm, p_ship_gwm, mdl.Coefficients.Estimate(1), mdl.Coefficients.Estimate(2), ...
    'VariableNames',{'r','p','Intercept','Slope'});
writetable(tbl_scatter1, scatter_xlsx, 'Sheet','Scatter_ShipleyVocab_GWMHFN');
writetable(tbl_meta1, scatter_xlsx, 'Sheet','Scatter_ShipleyVocab_GWMHFN','Range','E1');

% Figure 2: EstIQ_Shipley_Int_Bin (GWM-HFN)
selected_row = filtered_results_GWM(2);
X = selected_row.XS;
T = readtable(fullfile('../../Data/BGSP/','Behav_Data',[target_names{1} '.csv']),'Delimiter',',');
[isMatch, matchIdx] = ismember(ID, T.Match_ID);
Y = T{matchIdx(isMatch),2};
figure; out = scatplot(X, Y, [], [], [], [], [], 22); %#ok<NASGU>
colorbar('off'); hold on;
mdl = fitlm(X, Y);
x_range = linspace(min(X), max(X), 100);
y_fit = predict(mdl, x_range');
plot(x_range, y_fit, 'LineWidth', 2.5, 'Color', 'red');
[r_est_gwm, p_est_gwm] = corr(Y, X, 'type','Pearson');
text(0.05, 0.95, sprintf('r = %.3f', r_est_gwm), 'Units','normalized','FontSize',20,'Color','black','FontWeight','bold');
text(0.05, 0.90, sprintf('p = %.3g', p_est_gwm), 'Units','normalized','FontSize',20,'Color','black','FontWeight','bold');
title('Correlation of PLS1 score and  Estimated IQ (Shipley)', ...
    'FontSize',23,'FontWeight','bold','Units','normalized','Position',[0.5, 1.02]);
xlabel('PLS1 score (GWM-HFN)','FontSize',21,'FontWeight','bold');
ylabel('Shipley Vocabulary Task Score','FontSize',21,'FontWeight','bold');
xlim([-0.15, 0.11]); ylim([80, 140]);
xticks(-0.15:0.05:0.15); yticks(80:20:140);
set(gca,'FontSize',16,'FontWeight','bold'); grid on;
set(gcf,'Position',[1200, 600, 700, 600]);
print(fullfile('./Figures/','Corr_PLS1_EstIQ_Shiplep_GWM-HFN.png'),'-dpng','-r300'); hold off;

% === NEW EXPORT === Scatter source data (Sheet 2)
tbl_scatter2 = table(ID(isMatch), X, Y, 'VariableNames',{'SubjectID','PLS1_Score','EstIQ_Shipley_Int_Bin'});
tbl_meta2 = table(r_est_gwm, p_est_gwm, mdl.Coefficients.Estimate(1), mdl.Coefficients.Estimate(2), ...
    'VariableNames',{'r','p','Intercept','Slope'});
writetable(tbl_scatter2, scatter_xlsx, 'Sheet','Scatter_EstIQShipley_GWMHFN');
writetable(tbl_meta2, scatter_xlsx, 'Sheet','Scatter_EstIQShipley_GWMHFN','Range','E1');

% Figure 3: GM-GM Shipley_Vocab_Raw
load('Edge_BH_PLS_GMGM_BGSP.mat');  % results_GMGM
S = load('../../Data/BGSP/subjectIDs.mat'); ID = S.modified_subjects;
all_names = {results_GMGM.BH_name};
sel_mask = ismember(all_names, target_names);
filtered_results_GMGM = results_GMGM(sel_mask);
selected_row = filtered_results_GMGM(4);
X = selected_row.XS;
T = readtable(fullfile('../../Data/BGSP/','Behav_Data',[target_names{2} '.csv']),'Delimiter',',');
[isMatch, matchIdx] = ismember(ID, T.Match_ID);
Y = T{matchIdx(isMatch),2};
figure; out = scatplot(X, Y, [], [], [], [], [], 22); %#ok<NASGU>
colorbar('off'); hold on;
mdl = fitlm(X, Y);
x_range = linspace(min(X), max(X), 100);
y_fit = predict(mdl, x_range');
plot(x_range, y_fit, 'LineWidth', 2.5, 'Color', 'red');
[r_ship_gmgm, p_ship_gmgm] = corr(Y, X, 'type','Pearson');
text(0.05, 0.95, sprintf('r = %.3f', r_ship_gmgm), 'Units','normalized','FontSize',20,'Color','black','FontWeight','bold');
text(0.05, 0.90, sprintf('p = %.3g', p_ship_gmgm), 'Units','normalized','FontSize',20,'Color','black','FontWeight','bold');
title('Correlation of PLS1 score and Shipley vocabulary task', ...
    'FontSize',23,'FontWeight','bold','Units','normalized','Position',[0.5, 1.02]);
xlabel('PLS1 score (GM-GM)','FontSize',21,'FontWeight','bold');
ylabel('EstIQ (Shipley)','FontSize',21,'FontWeight','bold');
xlim([-0.15, 0.12]); ylim([20, 45]);
xticks(-0.15:0.05:0.12); yticks(20:5:45);
set(gca,'FontSize',16,'FontWeight','bold'); grid on;
set(gcf,'Position',[1200, 600, 700, 600]);
print(fullfile('./Figures/','Corr_PLS1_Shipley_Vocab_Raw_GM-GM.png'),'-dpng','-r300'); hold off;

% === NEW EXPORT === Scatter source data (Sheet 3)
tbl_scatter3 = table(ID(isMatch), X, Y, 'VariableNames',{'SubjectID','PLS1_Score','Shipley_Vocab_Raw'});
tbl_meta3 = table(r_ship_gmgm, p_ship_gmgm, mdl.Coefficients.Estimate(1), mdl.Coefficients.Estimate(2), ...
    'VariableNames',{'r','p','Intercept','Slope'});
writetable(tbl_scatter3, scatter_xlsx, 'Sheet','Scatter_ShipleyVocab_GMGM');
writetable(tbl_meta3, scatter_xlsx, 'Sheet','Scatter_ShipleyVocab_GMGM','Range','E1');

disp('GWM-HFN & GM-GM PLS scatterplots completed and source data exported.');

%% ========================================================================
%  PART D. Summary comparison (console)
% ========================================================================
fprintf('\n================= PART D: Summary Comparison =====================\n');
fprintf('%-3s %-35s | %-10s %-10s | %-10s %-10s\n', ...
    'Idx','Behavior','VarY_GWMHFN','p_perm','VarY_GMGM','p_perm');
fprintf(repmat('-',1,92)); fprintf('\n');
for i = BH_index
    fprintf('%-3d %-35s | %-10.4f %-10.4g | %-10.4f %-10.4g\n', ...
        i, results(i).BH_name, ...
        results(i).variance_Y, results(i).p, ...
        results_GMGM(i).variance_Y, results_GMGM(i).p);
end

%% ========================================================================
%  PART E. BrainNet export (two behaviors, conditional)
% ========================================================================
fprintf('\n================= PART E (Explicit BrainNet Export, No Loops) =================\n');

BSR_thr = 2;
z_thr   = 3;

BrainNetSummary = table();

% ------------------------------------------------------------------------
% 1) Behavior: EstIQ_Shipley_Int_Bin
%    Export: GWM-HFN ONLY
% ------------------------------------------------------------------------
beh_EstIQ = 'EstIQ_Shipley_Int_Bin';
idx_EstIQ_GWM = find(strcmp({results.BH_name}, beh_EstIQ), 1);

fprintf('\n--- Behavior: %s (GWM-HFN only) ---\n', beh_EstIQ);
w_gwm  = results(idx_EstIQ_GWM).weights(:);
BR_gwm = results(idx_EstIQ_GWM).boot_BR(:);
wz_gwm = zscore(w_gwm);
flag_gwm = (abs(wz_gwm) > z_thr) & (abs(BR_gwm) > BSR_thr);
edges_gwm = find(flag_gwm);

E_gwm = zeros(roi_num);
if ~isempty(edges_gwm)
    rg = r_all(edges_gwm); cg = c_all(edges_gwm);
    idx_lin = sub2ind([roi_num roi_num], rg, cg);
    E_gwm(idx_lin) = w_gwm(edges_gwm);
    E_gwm(sub2ind([roi_num roi_num], cg, rg)) = w_gwm(edges_gwm);
end

fn_edge_gwm = sprintf('%s_StableZgt%d_GWMHFN.edge', beh_EstIQ, z_thr);
dlmwrite(fn_edge_gwm, E_gwm, 'delimiter','\t','precision','%.6f');

deg = sum(E_gwm~=0,2);
node_out = node_base;
if any(deg), node_out(:,5) = MaxNodeSize * deg / max(deg); else, node_out(:,5)=0; end
fn_node_gwm = strrep(fn_edge_gwm,'.edge','.node');
fid = fopen(fn_node_gwm,'w');
for nn = 1:roi_num
    fprintf(fid,'%.6f\t%.6f\t%.6f\t%d\t%.6f\t%s\n', ...
        node_out(nn,1), node_out(nn,2), node_out(nn,3), ...
        node_out(nn,4), node_out(nn,5), node_labels{nn});
end
fclose(fid);

fprintf('  Exported (GWM-HFN): %s / %s (edges = %d)\n', fn_edge_gwm, fn_node_gwm, numel(edges_gwm));

% Append summary row
BrainNetSummary = [BrainNetSummary; ...
    table({beh_EstIQ},{'GWMHFN'},numel(edges_gwm),BSR_thr,z_thr, ...
    'VariableNames',{'Behavior','Type','EdgeCount','BSR_thr','Z_thr'})]; %#ok<AGROW>


% ------------------------------------------------------------------------
% 2) Behavior: Shipley_Vocab_Raw
%    Export: GWM-HFN + GM-GM + Consensus
% ------------------------------------------------------------------------
beh_Vocab = 'Shipley_Vocab_Raw';
idx_Vocab_GWM = find(strcmp({results.BH_name}, beh_Vocab), 1);
idx_Vocab_GMGM = find(strcmp({results_GMGM.BH_name}, beh_Vocab), 1);


fprintf('\n--- Behavior: %s (GWM-HFN + GM-GM + Consensus) ---\n', beh_Vocab);

% ---- GWM-HFN ----
w_gwm  = results(idx_Vocab_GWM).weights(:);
BR_gwm = results(idx_Vocab_GWM).boot_BR(:);
wz_gwm = zscore(w_gwm);
flag_gwm = (abs(wz_gwm) > z_thr) & (abs(BR_gwm) > BSR_thr);
edges_gwm = find(flag_gwm);

E_gwm = zeros(roi_num);
if ~isempty(edges_gwm)
    rg = r_all(edges_gwm); cg = c_all(edges_gwm);
    idx_lin = sub2ind([roi_num roi_num], rg, cg);
    E_gwm(idx_lin) = w_gwm(edges_gwm);
    E_gwm(sub2ind([roi_num roi_num], cg, rg)) = w_gwm(edges_gwm);
end

fn_edge_gwm = sprintf('%s_StableZgt%d_GWMHFN.edge', beh_Vocab, z_thr);
dlmwrite(fn_edge_gwm, E_gwm, 'delimiter','\t','precision','%.6f');

deg = sum(E_gwm~=0,2);
node_out = node_base;
if any(deg), node_out(:,5) = MaxNodeSize * deg / max(deg); else, node_out(:,5)=0; end
fn_node_gwm = strrep(fn_edge_gwm,'.edge','.node');
fid = fopen(fn_node_gwm,'w');
for nn = 1:roi_num
    fprintf(fid,'%.6f\t%.6f\t%.6f\t%d\t%.6f\t%s\n', ...
        node_out(nn,1), node_out(nn,2), node_out(nn,3), ...
        node_out(nn,4), node_out(nn,5), node_labels{nn});
end
fclose(fid);
fprintf('  Exported (GWM-HFN): %s / %s (edges = %d)\n', fn_edge_gwm, fn_node_gwm, numel(edges_gwm));
BrainNetSummary = [BrainNetSummary; ...
    table({beh_Vocab},{'GWMHFN'},numel(edges_gwm),BSR_thr,z_thr, ...
    'VariableNames',{'Behavior','Type','EdgeCount','BSR_thr','Z_thr'})]; %#ok<AGROW>

% ---- GM-GM ----

w_gmgm  = results_GMGM(idx_Vocab_GMGM).weights(:);
BR_gmgm = results_GMGM(idx_Vocab_GMGM).boot_BR(:);
wz_gmgm = zscore(w_gmgm);
flag_gmgm = (abs(wz_gmgm) > z_thr) & (abs(BR_gmgm) > BSR_thr);
edges_gmgm = find(flag_gmgm);

E_gmgm = zeros(roi_num);
if ~isempty(edges_gmgm)
    rg2 = r_all(edges_gmgm); cg2 = c_all(edges_gmgm);
    idx_lin2 = sub2ind([roi_num roi_num], rg2, cg2);
    E_gmgm(idx_lin2) = w_gmgm(edges_gmgm);
    E_gmgm(sub2ind([roi_num roi_num], cg2, rg2)) = w_gmgm(edges_gmgm);
end

fn_edge_gmgm = sprintf('%s_StableZgt%d_GMGM.edge', beh_Vocab, z_thr);
dlmwrite(fn_edge_gmgm, E_gmgm, 'delimiter','\t','precision','%.6f');

deg = sum(E_gmgm~=0,2);
node_out = node_base;
if any(deg), node_out(:,5) = MaxNodeSize * deg / max(deg); else, node_out(:,5)=0; end
fn_node_gmgm = strrep(fn_edge_gmgm,'.edge','.node');
fid = fopen(fn_node_gmgm,'w');
for nn = 1:roi_num
    fprintf(fid,'%.6f\t%.6f\t%.6f\t%d\t%.6f\t%s\n', ...
        node_out(nn,1), node_out(nn,2), node_out(nn,3), ...
        node_out(nn,4), node_out(nn,5), node_labels{nn});
end
fclose(fid);
fprintf('  Exported (GM-GM):   %s / %s (edges = %d)\n', fn_edge_gmgm, fn_node_gmgm, numel(edges_gmgm));
BrainNetSummary = [BrainNetSummary; ...
    table({beh_Vocab},{'GMGM'},numel(edges_gmgm),BSR_thr,z_thr, ...
    'VariableNames',{'Behavior','Type','EdgeCount','BSR_thr','Z_thr'})]; %#ok<AGROW>

% ---- Consensus (intersection) ----
if ~isempty(edges_gwm) && ~isempty(edges_gmgm)
    edges_intersect = intersect(edges_gwm, edges_gmgm);
else
    edges_intersect = [];
end

E_consensus = zeros(roi_num);
if ~isempty(edges_intersect)
    ri = r_all(edges_intersect); ci = c_all(edges_intersect);
    idx_linC = sub2ind([roi_num roi_num], ri, ci);
    mean_w = (w_gwm(edges_intersect) + w_gmgm(edges_intersect))/2;
    E_consensus(idx_linC) = mean_w;
    E_consensus(sub2ind([roi_num roi_num], ci, ri)) = mean_w;
end

fn_edge_cons = sprintf('%s_StableZgt%d_Consensus.edge', beh_Vocab, z_thr);
dlmwrite(fn_edge_cons, E_consensus, 'delimiter','\t','precision','%.6f');

deg = sum(E_consensus~=0,2);
node_out = node_base;
if any(deg), node_out(:,5) = MaxNodeSize * deg / max(deg); else, node_out(:,5)=0; end
fn_node_cons = strrep(fn_edge_cons,'.edge','.node');
fid = fopen(fn_node_cons,'w');
for nn = 1:roi_num
    fprintf(fid,'%.6f\t%.6f\t%.6f\t%d\t%.6f\t%s\n', ...
        node_out(nn,1), node_out(nn,2), node_out(nn,3), ...
        node_out(nn,4), node_out(nn,5), node_labels{nn});
end
fclose(fid);
fprintf('  Exported (Consensus): %s / %s (edges = %d)\n', fn_edge_cons, fn_node_cons, numel(edges_intersect));

BrainNetSummary = [BrainNetSummary; ...
    table({beh_Vocab},{'Consensus'},numel(edges_intersect),BSR_thr,z_thr, ...
    'VariableNames',{'Behavior','Type','EdgeCount','BSR_thr','Z_thr'})]; %#ok<AGROW>


% ------------------------------------------------------------------------
% Write summary table
% ------------------------------------------------------------------------
if ~isempty(BrainNetSummary)
    writetable(BrainNetSummary,'BrainNet_Export_Summary.csv');
    fprintf('\nBrainNet export summary saved: BrainNet_Export_Summary.csv\n');
else
    fprintf('\nNo BrainNet exports generated (summary empty).\n');
end

fprintf('\n================= END PART E (Explicit Export Completed) =====================\n');
%% ========================================================================
%  END
% ========================================================================
fprintf('\nAll analyses completed. Conditional BrainNet exports + comprehensive data exports generated.\n');