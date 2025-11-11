% ================================================================
% Age-Association Analysis (Schaefer100) for SALD:
%   - GWM-HFN (GM-WM-GM) vs GM-GM comparative modeling
%   - Linear vs Quadratic (AIC selection) per edge
%   - BH-FDR, network-level within/between distribution
%   - Overlap and effect-direction comparison between methods
%   - Peak age estimation for quadratic edges
%   - Optional GAM (global metric only)
%
% Assumptions:
%   * Subject ordering identical across matrices and demographics
%   * No missing subjects or reordering required
%   * Demographic .mat contains a table or struct with Age / Sex / Handedness / meanFD
%
% Author: (Your Name)
% Date: 2025-09-13
% ================================================================

clear; clc;

%% ===================== User Paths =====================
root_dir     = 'E:\Neuroimage\MyProject\GMWM_Network\Data\SALD';
file_gwmhfn  = fullfile(root_dir,'GWM_HFN_time1_sald_schaefer100.mat');
file_gmgm    = fullfile(root_dir,'GMGM_time1_sald_schaefer100.mat');
demo_csv     = fullfile(root_dir,'Anlysis_DemographicInfo_Schaefer.csv'); 
file_netidx  = fullfile(root_dir,'Schaefer100_7networkIndex.csv'); % length=100

out_dir      = fullfile(root_dir,'AgeAssociation_Schaefer100');
fig_dir      = fullfile(out_dir,'fig');
if ~exist(out_dir,'dir'); mkdir(out_dir); end
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

alpha_fdr = 0.05;
do_parallel = true;  % 是否使用 parfor 加速逐边拟合
save_edge_tables = true;

%% ===================== Load Connectivity Matrices =====================
S1 = load(file_gwmhfn);
S2 = load(file_gmgm);

% 自动识别 100x100xN 变量
fn1 = fieldnames(S1);
fn2 = fieldnames(S2);

mat_gwm = []; var_gwm = '';
for i=1:numel(fn1)
    v = S1.(fn1{i});
    if isnumeric(v) && ndims(v)==3 && size(v,1)==100 && size(v,2)==100
        mat_gwm = v; var_gwm = fn1{i}; break;
    end
end
if isempty(mat_gwm)
    error('未在 %s 中找到 100x100xN GWM-HFN 矩阵变量', file_gwmhfn);
end

mat_gmgm = []; var_gmgm = '';
for i=1:numel(fn2)
    v = S2.(fn2{i});
    if isnumeric(v) && ndims(v)==3 && size(v,1)==100 && size(v,2)==100
        mat_gmgm = v; var_gmgm = fn2{i}; break;
    end
end
if isempty(mat_gmgm)
    error('未在 %s 中找到 100x100xN GM-GM 矩阵变量', file_gmgm);
end

[n1,n2,N] = size(mat_gwm);
if any([n1 n2]~=100) || any(size(mat_gmgm,1:2)~=100) || size(mat_gmgm,3)~=N
    error('两个矩阵尺寸或被试数不一致。');
end
fprintf('Loaded: GWM-HFN var=%s, GM-GM var=%s, N=%d subjects\n', var_gwm, var_gmgm, N);

%% ===================== Load Demographics =====================
D = readtable(demo_csv,'TextType','string');
col_lower = lower(string(D.Properties.VariableNames));

col_alias = struct( ...
    'subid',  ["subid","subject","id","sub_id","participant_id"], ...
    'sex',    ["sex","gender"], ...
    'age',    ["age","years"], ...
    'hand',   ["hand","handedness"], ...
    'meanfd', ["meanfd","fd","fd_mean","mean_fd","mean_fd_jenkinson"] ...
);
get_col = @(cset) D.( D.Properties.VariableNames{ find( ismember(col_lower,cset),1,'first') } );

% 抓列（若缺失 Hand 列则默认全 Right）
try    SubID_raw = string(get_col(col_alias.subid)); catch, SubID_raw = strings(height(D),1); end
Age       = double(get_col(col_alias.age));
Sex_raw   = upper(string(get_col(col_alias.sex)));
try
    Hand_raw = upper(string(get_col(col_alias.hand)));
catch
    Hand_raw = repmat("RIGHT", height(D),1);
end
meanFD    = double(get_col(col_alias.meanfd));

if height(D) ~= N
    warning('人口学行数(%d)与连接矩阵被试数(%d)不一致，将截断对齐。', height(D), N);
    n_use = min(height(D), N);
    Age     = Age(1:n_use);
    Sex_raw = Sex_raw(1:n_use);
    Hand_raw= Hand_raw(1:n_use);
    meanFD  = meanFD(1:n_use);
    SubID_raw = SubID_raw(1:n_use);
    N = n_use;
end

SexM  = double(ismember(Sex_raw, ["M","MALE","1"]));
HandR = double(ismember(Hand_raw, ["R","RIGHT","1"]));
AgeZ  = (Age - mean(Age))/std(Age);
AgeZ2 = AgeZ.^2;
FDZ   = (meanFD - mean(meanFD))/std(meanFD);

fprintf('Demographics: N=%d, Age %.2f±%.2f, meanFD %.4f±%.4f, Male%%=%.1f, RightHand%%=%.1f\n',...
    N, mean(Age), std(Age), mean(meanFD), std(meanFD), 100*mean(SexM), 100*mean(HandR));
%% ===================== Design Matrices =====================
X_lin  = [AgeZ      SexM HandR FDZ];
X_quad = [AgeZ AgeZ2 SexM HandR FDZ];

varNames_lin  = {'AgeZ','SexM','HandR','FDZ','Y'};
varNames_quad = {'AgeZ','AgeZ2','SexM','HandR','FDZ','Y'};

%% ===================== Global Mean (Lower Triangle) =====================
maskLT = tril(true(100), -1);
idxLT  = find(maskLT);
Global_GWM = zeros(N,1);
Global_GM  = zeros(N,1);
for s=1:N
    M1 = mat_gwm(:,:,s);
    M2 = mat_gmgm(:,:,s);
    Global_GWM(s) = mean(M1(idxLT));
    Global_GM(s)  = mean(M2(idxLT));
end

% 线性 vs 二次
mdl_gwm_L = fitlm(X_lin,  Global_GWM, 'linear', 'VarNames',[varNames_lin(1:4) {'Global_GWM'}]);
mdl_gwm_Q = fitlm(X_quad, Global_GWM, 'linear', 'VarNames',[varNames_quad(1:5) {'Global_GWM'}]);
mdl_gm_L  = fitlm(X_lin,  Global_GM,  'linear', 'VarNames',[varNames_lin(1:4) {'Global_GM'}]);
mdl_gm_Q  = fitlm(X_quad, Global_GM,  'linear', 'VarNames',[varNames_quad(1:5) {'Global_GM'}]);

AIC_gwm_L = mdl_gwm_L.ModelCriterion.AIC; AIC_gwm_Q = mdl_gwm_Q.ModelCriterion.AIC;
AIC_gm_L  = mdl_gm_L.ModelCriterion.AIC; AIC_gm_Q  = mdl_gm_Q.ModelCriterion.AIC;

best_gwm = ternary(AIC_gwm_Q < AIC_gwm_L, 'Quadratic','Linear');
best_gm  = ternary(AIC_gm_Q  < AIC_gm_L,  'Quadratic','Linear');

fprintf('\nGlobal Mean Results:\n');
fprintf('  GWM-HFN: AIC(L)=%.2f  AIC(Q)=%.2f  Best=%s\n', AIC_gwm_L, AIC_gwm_Q, best_gwm);
fprintf('  GM-GM  : AIC(L)=%.2f  AIC(Q)=%.2f  Best=%s\n', AIC_gm_L,  AIC_gm_Q,  best_gm);

%% (Optional) GAM for Global GWM-HFN
try
    mdl_gam = fitrgam(table(Age,SexM,HandR,meanFD,Global_GWM), 'Global_GWM ~ s(Age) + SexM + HandR + meanFD');
    fprintf('GAM (Global GWM-HFN) trained. (仅探索，不参与后续统计)\n');
catch
    mdl_gam = [];
    warning('fitrgam 不可用或版本不支持，跳过 GAM。');
end

%% ===================== Edgewise Modeling (GWM-HFN) =====================
fprintf('\nEdgewise modeling: GWM-HFN (%d edges)...\n', numel(idxLT));
beta_age_gwm    = nan(100);   p_age_gwm    = nan(100); t_age_gwm    = nan(100);
beta_age2_gwm   = nan(100);   p_age2_gwm   = nan(100); t_age2_gwm   = nan(100);
model_used_gwm  = zeros(100); % 1=linear, 2=quadratic
AIC_L_gwm       = nan(100);   AIC_Q_gwm    = nan(100);
peak_age_gwm    = nan(100);   % raw years

[idx_i, idx_j] = find(maskLT);


for k=1:numel(idx_i)
    i = idx_i(k); j = idx_j(k);
    Y = squeeze(mat_gwm(i,j,:));
    mdlL = fitlm(X_lin,  Y, 'linear', 'VarNames', varNames_lin);
    mdlQ = fitlm(X_quad, Y, 'linear', 'VarNames', varNames_quad);
    AICL = mdlL.ModelCriterion.AIC;
    AICQ = mdlQ.ModelCriterion.AIC;
    if AICQ < AICL
        model_used_gwm(i,j)=2;
        AIC_Q_gwm(i,j)=AICQ; AIC_L_gwm(i,j)=AICL;
        c1 = mdlQ.Coefficients.Estimate('AgeZ');  p1 = mdlQ.Coefficients.pValue('AgeZ');  t1 = mdlQ.Coefficients.tStat('AgeZ');
        c2 = mdlQ.Coefficients.Estimate('AgeZ2'); p2 = mdlQ.Coefficients.pValue('AgeZ2'); t2 = mdlQ.Coefficients.tStat('AgeZ2');
        beta_age_gwm(i,j)=c1; p_age_gwm(i,j)=p1; t_age_gwm(i,j)=t1;
        beta_age2_gwm(i,j)=c2; p_age2_gwm(i,j)=p2; t_age2_gwm(i,j)=t2;
        if ~isnan(c2) && c2~=0
            peak_age_gwm(i,j) = -c1/(2*c2)*std(Age) + mean(Age);
        end
    else
        model_used_gwm(i,j)=1;
        AIC_Q_gwm(i,j)=AICQ; AIC_L_gwm(i,j)=AICL;
        c1 = mdlL.Coefficients.Estimate('AgeZ'); p1 = mdlL.Coefficients.pValue('AgeZ'); t1 = mdlL.Coefficients.tStat('AgeZ');
        beta_age_gwm(i,j)=c1; p_age_gwm(i,j)=p1; t_age_gwm(i,j)=t1;
    end
end


% 组合 p
p_comb_gwm = nan(100);
mask_lin_gwm  = (model_used_gwm==1);
mask_quad_gwm = (model_used_gwm==2);
p_comb_gwm(mask_lin_gwm)  = p_age_gwm(mask_lin_gwm);
p_comb_gwm(mask_quad_gwm) = p_age2_gwm(mask_quad_gwm);

% FDR
p_vec = p_comb_gwm(maskLT);
valid = ~isnan(p_vec);
m = sum(valid);
q_vec = nan(size(p_vec));
if m>0
    pv = p_vec(valid);
    [ps,ord] = sort(pv);
    ranks = (1:m)';
    q_sort = (m./ranks).*ps;
    q_sort = flipud(cummin(flipud(q_sort)));
    q_sort = min(q_sort,1);
    inv = zeros(m,1); inv(ord)=1:m;
    q_vec(valid) = q_sort(inv);
end
q_mat_gwm = nan(100); q_mat_gwm(maskLT)=q_vec;

sig_mask_gwm = (q_mat_gwm<alpha_fdr);
n_sig_gwm = nnz(sig_mask_gwm(maskLT));
fprintf('GWM-HFN significant edges (q<%.2f): %d / %d\n', alpha_fdr, n_sig_gwm, numel(idxLT));

% 系数矩阵（仅显著的保留）
matrix_age_gwm  = zeros(100);
matrix_age2_gwm = zeros(100);
sel_lin  = find(sig_mask_gwm & mask_lin_gwm & maskLT);
sel_quad = find(sig_mask_gwm & mask_quad_gwm & maskLT);
for k=1:numel(sel_lin)
    [i,j]=ind2sub([100 100], sel_lin(k));
    matrix_age_gwm(i,j)=beta_age_gwm(i,j);
end
for k=1:numel(sel_quad)
    [i,j]=ind2sub([100 100], sel_quad(k));
    matrix_age2_gwm(i,j)=beta_age2_gwm(i,j);
end

%% ===================== Edgewise Modeling (GM-GM) =====================
fprintf('\nEdgewise modeling: GM-GM (%d edges)...\n', numel(idxLT));
beta_age_gm    = nan(100);  p_age_gm    = nan(100); t_age_gm    = nan(100);
beta_age2_gm   = nan(100);  p_age2_gm   = nan(100); t_age2_gm   = nan(100);
model_used_gm  = zeros(100);
AIC_L_gm       = nan(100);  AIC_Q_gm    = nan(100);
peak_age_gm    = nan(100);

for k=1:numel(idx_i)
    i=idx_i(k); j=idx_j(k);
    Y = squeeze(mat_gmgm(i,j,:));
    mdlL = fitlm(X_lin, Y,'linear','VarNames',varNames_lin);
    mdlQ = fitlm(X_quad,Y,'linear','VarNames',varNames_quad);
    AICL=mdlL.ModelCriterion.AIC; AICQ=mdlQ.ModelCriterion.AIC;
    if AICQ<AICL
        model_used_gm(i,j)=2;
        AIC_Q_gm(i,j)=AICQ; AIC_L_gm(i,j)=AICL;
        c1=mdlQ.Coefficients.Estimate('AgeZ'); p1=mdlQ.Coefficients.pValue('AgeZ'); t1=mdlQ.Coefficients.tStat('AgeZ');
        c2=mdlQ.Coefficients.Estimate('AgeZ2');p2=mdlQ.Coefficients.pValue('AgeZ2'); t2=mdlQ.Coefficients.tStat('AgeZ2');
        beta_age_gm(i,j)=c1; p_age_gm(i,j)=p1; t_age_gm(i,j)=t1;
        beta_age2_gm(i,j)=c2; p_age2_gm(i,j)=p2; t_age2_gm(i,j)=t2;
        if ~isnan(c2) && c2~=0
            peak_age_gm(i,j) = -c1/(2*c2)*std(Age)+mean(Age);
        end
    else
        model_used_gm(i,j)=1;
        AIC_Q_gm(i,j)=AICQ; AIC_L_gm(i,j)=AICL;
        c1=mdlL.Coefficients.Estimate('AgeZ'); p1=mdlL.Coefficients.pValue('AgeZ'); t1=mdlL.Coefficients.tStat('AgeZ');
        beta_age_gm(i,j)=c1; p_age_gm(i,j)=p1; t_age_gm(i,j)=t1;
    end
end
  

% 组合 p 与 FDR
p_comb_gm = nan(100);
p_comb_gm(model_used_gm==1) = p_age_gm(model_used_gm==1);
p_comb_gm(model_used_gm==2) = p_age2_gm(model_used_gm==2);

p_vec2 = p_comb_gm(maskLT);
valid2 = ~isnan(p_vec2);
m2 = sum(valid2);
q_vec2 = nan(size(p_vec2));
if m2>0
    pv2 = p_vec2(valid2);
    [ps2,ord2] = sort(pv2);
    ranks2 = (1:m2)';
    q_sort2 = (m2./ranks2).*ps2;
    q_sort2 = flipud(cummin(flipud(q_sort2)));
    q_sort2 = min(q_sort2,1);
    inv2 = zeros(m2,1); inv2(ord2)=1:m2;
    q_vec2(valid2)=q_sort2(inv2);
end
q_mat_gm = nan(100); q_mat_gm(maskLT)=q_vec2;
sig_mask_gm = q_mat_gm<alpha_fdr;
n_sig_gm = nnz(sig_mask_gm(maskLT));
fprintf('GM-GM significant edges (q<%.2f): %d / %d\n', alpha_fdr, n_sig_gm, numel(idxLT));

matrix_age_gm  = zeros(100);
matrix_age2_gm = zeros(100);
sel_lin_gm  = find(sig_mask_gm & model_used_gm==1 & maskLT);
sel_quad_gm = find(sig_mask_gm & model_used_gm==2 & maskLT);
for k=1:numel(sel_lin_gm)
    [i,j]=ind2sub([100 100], sel_lin_gm(k));
    matrix_age_gm(i,j)=beta_age_gm(i,j);
end
for k=1:numel(sel_quad_gm)
    [i,j]=ind2sub([100 100], sel_quad_gm(k));
    matrix_age2_gm(i,j)=beta_age2_gm(i,j);
end

%% ===================== Overlap & Direction Consistency =====================
lin_mask_gwm_sig = (matrix_age_gwm~=0);
lin_mask_gm_sig  = (matrix_age_gm~=0);
quad_mask_gwm_sig= (matrix_age2_gwm~=0);
quad_mask_gm_sig = (matrix_age2_gm~=0);

% 合并（线性+二次）
all_sig_gwm = (lin_mask_gwm_sig | quad_mask_gwm_sig);
all_sig_gm  = (lin_mask_gm_sig  | quad_mask_gm_sig);

overlap_all = all_sig_gwm & all_sig_gm;
n_overlap   = nnz(overlap_all(maskLT));
jaccard_all = n_overlap / ( nnz(all_sig_gwm(maskLT)) + nnz(all_sig_gm(maskLT)) - n_overlap );

% 方向一致性（仅对两方法都显著的线性项：若某边在一侧是线性在另一侧是二次→用对应项的"AgeZ"系数；二次模型中 AgeZ 为一次项系数）
% 提取 AgeZ 系数矩阵：线性已经在 beta_age_*；二次也存储在 beta_age_*中
coef_gwm_primary = beta_age_gwm;  % AgeZ 系数（对所有边模型：线性/二次都有 AgeZ）
coef_gm_primary  = beta_age_gm;

overlap_idx = find(overlap_all & maskLT);
agree_count = 0; opposite_count=0;
for u=1:numel(overlap_idx)
    [i,j] = ind2sub([100 100], overlap_idx(u));
    c1 = coef_gwm_primary(i,j);
    c2 = coef_gm_primary(i,j);
    if c1*c2 > 0
        agree_count = agree_count+1;
    elseif c1*c2 < 0
        opposite_count = opposite_count+1;
    end
end
direction_agreement_rate = agree_count / max(1, numel(overlap_idx));

fprintf('\nOverlap Summary (All Significant Edges):\n');
fprintf('  GWM-HFN sig edges: %d\n', nnz(all_sig_gwm(maskLT)));
fprintf('  GM-GM   sig edges: %d\n', nnz(all_sig_gm(maskLT)));
fprintf('  Overlap          : %d\n', n_overlap);
fprintf('  Jaccard          : %.3f\n', jaccard_all);
fprintf('  Direction agreement (AgeZ sign) among overlaps: %.2f%% (%d/%d)\n', ...
    100*direction_agreement_rate, agree_count, max(1,numel(overlap_idx)));

%% ===================== Network-Level (Schaefer 7) =====================
net_index = readmatrix(file_netidx); % length=100
if numel(net_index)~=100
    error('网络索引长度不是100。');
end
[~,~,net_index_remap] = unique(net_index);
K = numel(unique(net_index_remap));

% combined significant edges（GWM-HFN）
sig_comb_gwm_LT = all_sig_gwm & maskLT;
within_mask = false(100);
for k=1:K
    nodes = find(net_index_remap==k);
    tmp = false(100);
    tmp(nodes,nodes)=true;
    within_mask = within_mask | (tmp & maskLT);
end
between_mask = maskLT & ~within_mask;

Num_within_gwm  = nnz(sig_comb_gwm_LT & within_mask);
Num_between_gwm = nnz(sig_comb_gwm_LT & between_mask);
Den_total = 100*99/2;
Den_within = 0;
for k=1:K
    nk = sum(net_index_remap==k);
    Den_within = Den_within + nk*(nk-1)/2;
end
Den_between = Den_total - Den_within;

Prop_within_gwm  = Num_within_gwm / Den_within;
Prop_between_gwm = Num_between_gwm / Den_between;

fprintf('\nNetwork Distribution (GWM-HFN combined):\n');
fprintf('  Within: %d / %d = %.4f\n', Num_within_gwm,  Den_within,  Prop_within_gwm);
fprintf('  Between:%d / %d = %.4f\n', Num_between_gwm, Den_between, Prop_between_gwm);

%% ===================== Save Edge Tables =====================
if save_edge_tables
    % 生成边列表 (GWM-HFN)
    Node1 = idx_i; Node2=idx_j;
    ModelUsed = model_used_gwm(maskLT);
    Beta_Age  = beta_age_gwm(maskLT);
    Beta_Age2 = beta_age2_gwm(maskLT);
    P_Age     = p_age_gwm(maskLT);
    P_Age2    = p_age2_gwm(maskLT);
    Q_value   = q_mat_gwm(maskLT);
    PeakAge   = peak_age_gwm(maskLT);
    T_edge_gwm = table(Node1,Node2,ModelUsed,Beta_Age,Beta_Age2,P_Age,P_Age2,Q_value,PeakAge);
    writetable(T_edge_gwm, fullfile(out_dir,'Edges_GWMHFN_Schaefer100_Age.csv'));

    % GM-GM
    ModelUsed_gm = model_used_gm(maskLT);
    Beta_Age_gm  = beta_age_gm(maskLT);
    Beta_Age2_gm = beta_age2_gm(maskLT);
    P_Age_gm     = p_age_gm(maskLT);
    P_Age2_gm    = p_age2_gm(maskLT);
    Q_value_gm   = q_mat_gm(maskLT);
    PeakAge_gm   = peak_age_gm(maskLT);
    T_edge_gm = table(Node1,Node2,ModelUsed_gm,Beta_Age_gm,Beta_Age2_gm,P_Age_gm,P_Age2_gm,Q_value_gm,PeakAge_gm);
    writetable(T_edge_gm, fullfile(out_dir,'Edges_GMGM_Schaefer100_Age.csv'));

    % Overlap summary
    OverlapSummary = table(...
        nnz(all_sig_gwm(maskLT)), nnz(all_sig_gm(maskLT)), n_overlap, jaccard_all, ...
        direction_agreement_rate, ...
        'VariableNames', {'Sig_GWMHFN','Sig_GMGM','Overlap','Jaccard','DirAgreement'});
    writetable(OverlapSummary, fullfile(out_dir,'Overlap_GWMHFN_vs_GMGM.csv'));
end

%% ===================== Figures (Quick Matrix Heatmaps) =====================
f1 = figure('Color','w'); imagesc(matrix_age_gwm); axis square; colorbar;
title('GWM-HFN Linear AgeZ Coeff (q<0.05)');
print(f1, fullfile(fig_dir,'GWMHFN_LinearAgeCoeff.tiff'),'-dtiff','-r300');

f2 = figure('Color','w'); imagesc(matrix_age2_gwm); axis square; colorbar;
title('GWM-HFN Quadratic AgeZ^2 Coeff (q<0.05)');
print(f2, fullfile(fig_dir,'GWMHFN_QuadraticAgeCoeff.tiff'),'-dtiff','-r300');

f3 = figure('Color','w'); imagesc(matrix_age_gm); axis square; colorbar;
title('GM-GM Linear AgeZ Coeff (q<0.05)');
print(f3, fullfile(fig_dir,'GMGM_LinearAgeCoeff.tiff'),'-dtiff','-r300');

f4 = figure('Color','w'); imagesc(matrix_age2_gm); axis square; colorbar;
title('GM-GM Quadratic AgeZ^2 Coeff (q<0.05)');
print(f4, fullfile(fig_dir,'GMGM_QuadraticAgeCoeff.tiff'),'-dtiff','-r300');

%% ===================== Save Structured Results =====================
Results = struct();
Results.Meta = struct('N',N,'alpha_fdr',alpha_fdr,'atlas','Schaefer100',...
    'file_gwmhfn',file_gwmhfn,'file_gmgm',file_gmgm,'demo',file_demo);
Results.Global = struct(...
    'GWM_L',mdl_gwm_L,'GWM_Q',mdl_gwm_Q,'Best_GWM',best_gwm,...
    'GM_L', mdl_gm_L, 'GM_Q', mdl_gm_Q, 'Best_GM', best_gm,...
    'GAM_Global_GWM', mdl_gam);
Results.GWMHFN = struct(...
    'beta_age',beta_age_gwm,'beta_age2',beta_age2_gwm,...
    'p_age',p_age_gwm,'p_age2',p_age2_gwm,'q_value',q_mat_gwm,...
    'model_used',model_used_gwm,'peak_age_years',peak_age_gwm,...
    'matrix_age_sig',matrix_age_gwm,'matrix_age2_sig',matrix_age2_gwm);
Results.GMGM = struct(...
    'beta_age',beta_age_gm,'beta_age2',beta_age2_gm,...
    'p_age',p_age_gm,'p_age2',p_age2_gm,'q_value',q_mat_gm,...
    'model_used',model_used_gm,'peak_age_years',peak_age_gm,...
    'matrix_age_sig',matrix_age_gm,'matrix_age2_sig',matrix_age2_gm);
Results.Compare = struct(...
    'all_sig_gwm',all_sig_gwm,'all_sig_gm',all_sig_gm,...
    'overlap_mask',overlap_all,'jaccard',jaccard_all,...
    'direction_agreement_rate',direction_agreement_rate);
Results.NetworkDist = struct(...
    'Within_Count',Num_within_gwm,'Between_Count',Num_between_gwm,...
    'Within_Den',Den_within,'Between_Den',Den_between,...
    'Within_Prop',Prop_within_gwm,'Between_Prop',Prop_between_gwm);

save(fullfile(out_dir,'AgeEffects_Schaefer100_GWMHFN_vs_GMGM.mat'),'Results','-v7.3');

fprintf('\nPipeline finished. Outputs saved to: %s\n', out_dir);

%% ===================== Helper =====================
function out = ternary(cond,a,b)
if cond; out=a; else; out=b; end
end