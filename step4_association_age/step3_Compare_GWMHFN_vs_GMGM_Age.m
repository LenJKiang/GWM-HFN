% step3_Compare_GWMHFN_vs_GMGM_Age.m
% Comprehensive comparison of age effects between:
%   (1) GWM-HFN (GM-WM-GM derived network)
%   (2) Classical GM-GM functional connectivity
%
% Analytic components:
%   A. Edge-level overlap and directional / shape concordance
%   B. Broad directional grouping (decline-like vs increase-like)
%   C. Detailed overlap table export
%   D. Category heatmap (significance class) using lc_netplot
%   E. Scatter comparisons for (i) negative linear overlaps and (ii) inverted-U quadratic overlaps
%   F. Network-level distribution of significant edges (within / between networks)
%   G. Subject-level method agreement and its age association
%   H. Consolidated outputs (.csv, .mat, figures)
%
% NOTE (per user request):
%   - Plot aesthetics (styling, file names, formats) are kept intact.
%   - Code structure, comments, and logic have been streamlined.
%   - All comments standardized to English.
%
% Requirements (must exist beforehand):
%   AgeEffects_GWMHFN_SALD.mat  -> struct AgeEffects (fields in Edgewise)
%   AgeEffects_GMGM_SALD.mat    -> struct AgeEffects (fields in Edgewise)
%   FCMData_GWMHFN.mat          -> raw 90x90xN matrix fcm_GWG_cstd_prod
%   FCMData_GMGM.mat            -> raw 90x90xN matrix Z_GM
%   AAL_7networkIndex.csv       -> 90×1 network index (1..7)
%
% Author: Zeqiang LINLI
% Date: 2025-08-14

clear; clc;

%% ===================== (0) User Configuration =====================
data_root        = 'E:\Neuroimage\MyProject\GMWM_Network\Data\SALD';
age_dir          = 'E:\Neuroimage\MyProject\GMWM_Network\Code_R1\step4_association_age';
network_idx_file = 'E:\Neuroimage\MyProject\GMWM_Network\Data\AAL_7networkIndex.csv';
fig_dir          = fullfile('E:\Neuroimage\MyProject\GMWM_Network\Code_R1\visualization','Figure6_Age');
out_dir          = fullfile(age_dir,'comparison_GWMHFN_vs_GMGM');

alpha_fdr        = 0.05;
do_save_figures  = true;   % Global flag for saving figures
save_fig         = true;   % (Used in original plotting section, keep for compatibility)
useSizeEncoding  = true;   % Placeholder flag (not used for scaling in current version)

% Create directories if missing
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
if ~exist(out_dir,'dir'); mkdir(out_dir); end

%% ===================== (1) Load Age Effect Structures =====================
file_gwm = fullfile(age_dir,'AgeEffects_GWMHFN_SALD.mat');
file_gm  = fullfile(age_dir,'AgeEffects_GMGM_SALD.mat');

assert(exist(file_gwm,'file')==2, 'Missing file: %s', file_gwm);
assert(exist(file_gm,'file')==2,  'Missing file: %s', file_gm);

Sgwm = load(file_gwm); Age_GWM = Sgwm.AgeEffects;
Sgm  = load(file_gm);  Age_GM  = Sgm.AgeEffects;

% Required edgewise fields
requiredFields = {'q_value','model_used','beta_age','beta_age2','t_age','t_age2'};
assert(all(ismember(requiredFields, fieldnames(Age_GWM.Edgewise))), 'Missing fields in GWM-HFN AgeEffects.Edgewise');
assert(all(ismember(requiredFields, fieldnames(Age_GM.Edgewise))),  'Missing fields in GM-GM AgeEffects.Edgewise');

% Dimension sanity
assert(isequal(size(Age_GWM.Edgewise.q_value),[90 90]), 'GWM-HFN matrix must be 90x90');
assert(isequal(size(Age_GM.Edgewise.q_value),[90 90]),  'GM-GM matrix must be 90x90');

%% ===================== (2) Extract Core Edgewise Data =====================
% Abbreviations:
%   model_used: 1 = linear, 2 = quadratic
%   beta_age  : linear age effect (z-scored age)
%   beta_age2 : quadratic age effect (z-scored age^2)
%   t_age / t_age2 : t-statistics for linear / quadratic terms

q_gwm      = Age_GWM.Edgewise.q_value;
model_gwm  = Age_GWM.Edgewise.model_used;
betaL_gwm  = Age_GWM.Edgewise.beta_age;
betaQ_gwm  = Age_GWM.Edgewise.beta_age2;
t_lin_gwm  = Age_GWM.Edgewise.t_age;
t_quad_gwm = Age_GWM.Edgewise.t_age2;

q_gm       = Age_GM.Edgewise.q_value;
model_gm   = Age_GM.Edgewise.model_used;
betaL_gm   = Age_GM.Edgewise.beta_age;
betaQ_gm   = Age_GM.Edgewise.beta_age2;
t_lin_gm   = Age_GM.Edgewise.t_age;
t_quad_gm  = Age_GM.Edgewise.t_age2;

% Binary significance masks
sig_gwm = q_gwm < alpha_fdr;
sig_gm  = q_gm  < alpha_fdr;

% Lower triangle indices (unique undirected edges)
maskLT         = tril(true(90), -1);
[idx_i, idx_j] = find(maskLT);
n_edges        = numel(idx_i);

% Convenient lower-triangle masks
sig_gwm_LT = sig_gwm & maskLT;
sig_gm_LT  = sig_gm  & maskLT;

%% ===================== (3) Edge-Level Overlap Statistics =====================
% Overlap counts
n_sig_gwm      = nnz(sig_gwm_LT);
n_sig_gm       = nnz(sig_gm_LT);
overlap_sig_LT = sig_gwm_LT & sig_gm_LT;
n_overlap_sig  = nnz(overlap_sig_LT);

% Jaccard = |A ∩ B| / |A ∪ B|
Jaccard       = n_overlap_sig / (n_sig_gwm + n_sig_gm - n_overlap_sig);
Percent_gmgm  = n_overlap_sig / max(1,n_sig_gm);  % Overlap proportion relative to GM-GM significant

fprintf('\n=== (1) Edge-level significance overlap (q < %.2f) ===\n', alpha_fdr);
fprintf('GWM-HFN significant edges : %d\n', n_sig_gwm);
fprintf('GM-GM   significant edges : %d\n', n_sig_gm);
fprintf('Overlap significant edges : %d\n', n_overlap_sig);
fprintf('Jaccard similarity        : %.4f\n', Jaccard);
fprintf('Proportion overlap / GM-GM: %.4f\n', Percent_gmgm);
% GWM-HFN significant edges : 2044
% GM-GM   significant edges : 1550
% Overlap significant edges : 938
% Jaccard similarity        : 0.3532
% Proportion overlap / GM-GM: 0.6052
%% ===================== (4) Linear Direction Concordance (Overlap & Both Linear) =====================
lin_overlap_mask = overlap_sig_LT & (model_gwm==1) & (model_gm==1);
n_lin_overlap    = nnz(lin_overlap_mask);

lin_beta_gwm     = betaL_gwm(lin_overlap_mask);
lin_beta_gm      = betaL_gm(lin_overlap_mask);

same_direction   = sum( (lin_beta_gwm > 0 & lin_beta_gm > 0) | (lin_beta_gwm < 0 & lin_beta_gm < 0) );
opp_direction    = n_lin_overlap - same_direction;

fprintf('\nLinear overlap edges: %d\n', n_lin_overlap);
fprintf('  Same sign directions   : %d (%.1f%%)\n', same_direction, 100*same_direction/max(1,n_lin_overlap));
fprintf('  Opposite sign directions: %d (%.1f%%)\n', opp_direction,  100*opp_direction/max(1,n_lin_overlap));
% Linear overlap edges: 412
%   Same sign directions   : 411 (99.8%)
%   Opposite sign directions: 1 (0.2%)
%% ===================== (5) Quadratic Shape Concordance (Overlap & Both Quadratic) =====================
quad_overlap_mask = overlap_sig_LT & (model_gwm==2) & (model_gm==2);
n_quad_overlap    = nnz(quad_overlap_mask);

quad_beta_gwm     = betaQ_gwm(quad_overlap_mask);
quad_beta_gm      = betaQ_gm(quad_overlap_mask);

% Shape: >0 = U-shape; <0 = inverted-U
isU_gwm    = quad_beta_gwm > 0;
isU_gm     = quad_beta_gm  > 0;
same_shape = sum( (isU_gwm & isU_gm) | (~isU_gwm & ~isU_gm) );
diff_shape = n_quad_overlap - same_shape;

fprintf('\nQuadratic overlap edges: %d\n', n_quad_overlap);
fprintf('  Same shape (U/U or InvU/InvU): %d (%.1f%%)\n', same_shape, 100*same_shape/max(1,n_quad_overlap));
fprintf('  Different shape               : %d (%.1f%%)\n', diff_shape, 100*diff_shape/max(1,n_quad_overlap));
% Quadratic overlap edges: 287
%   Same shape (U/U or InvU/InvU): 287 (100.0%)
%   Different shape               : 0 (0.0%)
%% ===================== (6) Broad Directional Concordance =====================
% Group definitions:
%   Group A (Decline-like)  : Linear negative, Inverted-U (betaQ<0)
%   Group B (Increase-like) : Linear positive, U-shape    (betaQ>0)
% Classification per method only on significant edges.

beta_type_gwm = strings(90,90);
beta_type_gm  = strings(90,90);

% GWM-HFN classification
lin_mask_gwm_all  = (model_gwm==1) & sig_gwm;
quad_mask_gwm_all = (model_gwm==2) & sig_gwm;
beta_type_gwm(lin_mask_gwm_all  & betaL_gwm>0) = "LIN_POS";
beta_type_gwm(lin_mask_gwm_all  & betaL_gwm<0) = "LIN_NEG";
beta_type_gwm(quad_mask_gwm_all & betaQ_gwm>0) = "QUAD_U";
beta_type_gwm(quad_mask_gwm_all & betaQ_gwm<0) = "QUAD_INVU";

% GM-GM classification
lin_mask_gm_all  = (model_gm==1) & sig_gm;
quad_mask_gm_all = (model_gm==2) & sig_gm;
beta_type_gm(lin_mask_gm_all  & betaL_gm>0) = "LIN_POS";
beta_type_gm(lin_mask_gm_all  & betaL_gm<0) = "LIN_NEG";
beta_type_gm(quad_mask_gm_all & betaQ_gm>0) = "QUAD_U";
beta_type_gm(quad_mask_gm_all & betaQ_gm<0) = "QUAD_INVU";

GroupA = ["LIN_NEG","QUAD_INVU"];
GroupB = ["LIN_POS","QUAD_U"];

broad_overlap_mask  = overlap_sig_LT;
n_broad_overlap     = nnz(broad_overlap_mask);

tmpA = ismember(beta_type_gwm, GroupA) & ismember(beta_type_gm, GroupA);
tmpB = ismember(beta_type_gwm, GroupB) & ismember(beta_type_gm, GroupB);
broad_consistent_mask   = (tmpA | tmpB) & broad_overlap_mask;
n_broad_consistent      = nnz(broad_consistent_mask);
n_broad_inconsistent    = n_broad_overlap - n_broad_consistent;

fprintf('\nBroad directional overlap (all overlapping significant edges: %d)\n', n_broad_overlap);
fprintf('  Broadly consistent : %d (%.1f%%)\n', n_broad_consistent, 100*n_broad_consistent/max(1,n_broad_overlap));
fprintf('  Broadly inconsistent: %d (%.1f%%)\n', n_broad_inconsistent, 100*n_broad_inconsistent/max(1,n_broad_overlap));
% Broad directional overlap (all overlapping significant edges: 938)
%   Broadly consistent : 921 (98.2%)
%   Broadly inconsistent: 17 (1.8%)
%% ===================== (7) Assemble & Export Edge-Level Summary =====================
EdgeOverlapSummary = table( ...
    n_sig_gwm, n_sig_gm, n_overlap_sig, Jaccard, ...
    n_lin_overlap, same_direction, opp_direction, ...
    n_quad_overlap, same_shape, diff_shape, ...
    n_broad_overlap, n_broad_consistent, n_broad_inconsistent, ...
    'VariableNames', {'Sig_GWMHFN','Sig_GMGM','Sig_Overlap','Jaccard', ...
    'Linear_Overlap','Linear_SameDir','Linear_OppDir', ...
    'Quadratic_Overlap','Quadratic_SameShape','Quadratic_DiffShape', ...
    'Broad_Overlap','Broad_Consistent','Broad_Inconsistent'});

% Detailed overlapping edges
overlap_idx = find(overlap_sig_LT);
[ov_r, ov_c] = ind2sub([90,90], overlap_idx);
class_gwm = beta_type_gwm(overlap_idx);
class_gm  = beta_type_gm(overlap_idx);

betaL_gwm_full = betaL_gwm(overlap_idx);
betaL_gm_full  = betaL_gm(overlap_idx);
betaQ_gwm_full = betaQ_gwm(overlap_idx);
betaQ_gm_full  = betaQ_gm(overlap_idx);
broad_flag_vec = broad_consistent_mask(overlap_idx);

DetailedOverlap = table(ov_r, ov_c, ...
    class_gwm(:), class_gm(:), ...
    betaL_gwm_full(:), betaL_gm_full(:), betaQ_gwm_full(:), betaQ_gm_full(:), ...
    broad_flag_vec(:), ...
    'VariableNames', {'Node1','Node2','Class_GWMHFN','Class_GMGM', ...
    'BetaL_GWMHFN','BetaL_GMGM','BetaQ_GWMHFN','BetaQ_GMGM','BroadConsistent'});

writetable(EdgeOverlapSummary, fullfile(out_dir,'EdgeOverlapSummary.csv'));
writetable(DetailedOverlap,     fullfile(out_dir,'DetailedOverlapEdges.csv'));

%% ===================== (8) Significance Category Heatmap (lc_netplot) =====================
% Categories:
%   1 = Neither
%   2 = Only GM-GM
%   3 = Only GWM-HFN
%   4 = Both

sig_gwm = q_gwm < alpha_fdr;
sig_gwm = sig_gwm + sig_gwm'
sig_gm  = q_gm  < alpha_fdr;
sig_gm = sig_gm + sig_gm'

overlap_mask = sig_gwm & sig_gm;
only_gwm     = sig_gwm & ~sig_gm;
only_gm      = ~sig_gwm & sig_gm;


Category = ones(90,90,'uint8');
Category(only_gm)      = 2;
Category(only_gwm)     = 3;
Category(overlap_mask) = 4;

netIndex = readmatrix(network_idx_file);
assert(numel(netIndex)==90, 'Network index must have 90 entries.');
legends = {'VN','SMN','AN','LB','FPN','DMN','BGN'};

figure('Position',[100 100 420 420],'Color','w');
lc_netplot('-n', Category, '-ni', netIndex, '-il', 1, '-lg', legends);
axis square;

% Bold network labels (post-draw tweak)
drawnow;
ax   = gca;
txts = findall(ax,'Type','text');
legSet = string(legends);
for t = txts'
    if iscell(t.String), s = string(t.String{1}); else, s = string(t.String); end
    if ismember(strtrim(s), legSet)
        t.FontWeight = 'bold';
        t.FontSize   = 12;
    end
end

cmap = [0.95 0.95 0.95; [201,203,163]/255; [255,225,168]/255; [226,109,92]/255];
colormap(cmap); caxis([1 4]);
colorbar
title(sprintf('Significance Categories (P_{FDR}<%.2f)',alpha_fdr), ...
    'FontSize',13,'FontWeight','bold','Units','normalized','Position',[0.5,0.92,0]);

if save_fig
    print(fullfile(fig_dir,'Figure_OverlapCategoryHeatmap_lcnet.tiff'),'-dtiff','-r300');
    print(fullfile(fig_dir,'Figure_OverlapCategoryHeatmap_lcnet.svg'), '-dsvg');
end

%% ===================== (9) Scatter: Negative Linear & Inverted-U Quadratic Overlaps =====================
lin_neg_overlap = overlap_mask & (model_gwm==1) & (model_gm==1) & ...
                  (betaL_gwm<0) & (betaL_gm<0) & maskLT;
[idx_i_lin, idx_j_lin] = find(lin_neg_overlap);

quad_invU_overlap = overlap_mask & (model_gwm==2) & (model_gm==2) & ...
                    (betaQ_gwm<0) & (betaQ_gm<0) & maskLT;
[idx_i_q, idx_j_q] = find(quad_invU_overlap);

% Build vectors (negative linear)
if isempty(idx_i_lin)
    t_lin_vec_gwm = [];
    t_lin_vec_gm  = [];
else
    nL = numel(idx_i_lin);
    t_lin_vec_gwm = zeros(nL,1);
    t_lin_vec_gm  = zeros(nL,1);
    for k = 1:nL
        t_lin_vec_gwm(k) = t_lin_gwm(idx_i_lin(k), idx_j_lin(k));
        t_lin_vec_gm(k)  = t_lin_gm (idx_i_lin(k), idx_j_lin(k));
    end
end

% Build vectors (inverted-U quadratic)
if isempty(idx_i_q)
    t_quad_vec_gwm = [];
    t_quad_vec_gm  = [];
else
    nQ = numel(idx_i_q);
    t_quad_vec_gwm = zeros(nQ,1);
    t_quad_vec_gm  = zeros(nQ,1);
    for k = 1:nQ
        t_quad_vec_gwm(k) = t_quad_gwm(idx_i_q(k), idx_j_q(k));
        t_quad_vec_gm(k)  = t_quad_gm (idx_i_q(k), idx_j_q(k));
    end
end

% Safe correlations
if numel(t_lin_vec_gm) >= 2
    validMask = ~(isnan(t_lin_vec_gm) | isnan(t_lin_vec_gwm));
    if sum(validMask)>=2
        [r_lin_t,p_lin_t] = corr(t_lin_vec_gm(validMask), t_lin_vec_gwm(validMask),'type','Pearson');
    else
        r_lin_t = NaN; p_lin_t = NaN;
    end
else
    r_lin_t = NaN; p_lin_t = NaN;
end

if numel(t_quad_vec_gm) >= 2
    validMaskQ = ~(isnan(t_quad_vec_gm) | isnan(t_quad_vec_gwm));
    if sum(validMaskQ)>=2
        [r_quad_t,p_quad_t] = corr(t_quad_vec_gm(validMaskQ), t_quad_vec_gwm(validMaskQ),'type','Pearson');
    else
        r_quad_t = NaN; p_quad_t = NaN;
    end
else
    r_quad_t = NaN; p_quad_t = NaN;
end

% Fisher Z (store both z value and approximate p)
if isnan(r_lin_t) || numel(t_lin_vec_gm) < 4
    r_lin_z = NaN; p_lin_z = NaN;
else
    zf_lin = atanh(r_lin_t);
    z_stat_lin = zf_lin * sqrt(numel(t_lin_vec_gm)-3);
    p_lin_z = 2*(1-normcdf(abs(z_stat_lin),0,1));
    r_lin_z = zf_lin;
end
if isnan(r_quad_t) || numel(t_quad_vec_gm) < 4
    r_quad_z = NaN; p_quad_z = NaN;
else
    zf_quad = atanh(r_quad_t);
    z_stat_quad = zf_quad * sqrt(numel(t_quad_vec_gm)-3);
    p_quad_z = 2*(1-normcdf(abs(z_stat_quad),0,1));
    r_quad_z = zf_quad;
end

% Combined panel figure
figure('Color','w','Position',[100 100 1100 480]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Panel A: Negative linear
nexttile;
if isempty(t_lin_vec_gm)
    text(0.5,0.5,'No negative linear overlap edges','HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
    axis off
else
    scatter(t_lin_vec_gm, t_lin_vec_gwm, 40, 'MarkerFaceColor',[0.15 0.35 0.75], ...
        'MarkerEdgeColor','k','LineWidth',0.4,'MarkerFaceAlpha',0.75); hold on;
    mn = min([t_lin_vec_gm; t_lin_vec_gwm]); mx = max([t_lin_vec_gm; t_lin_vec_gwm]);
    if mx==mn, mn = mn - 0.5; mx = mx + 0.5; end
    pad = 0.05*(mx-mn);
    xlim([mn-pad mx+pad]); ylim([mn-pad mx+pad]);
    plot([mn mx],[mn mx],'k--','LineWidth',1.1);
    if numel(t_lin_vec_gm) >= 2
        mdlL = fitlm(t_lin_vec_gm, t_lin_vec_gwm);
        xf = linspace(xlim(gca)*[1;0], xlim(gca)*[0;1], 200);
        y_fit = predict(mdlL, xf');
        plot(xf, y_fit, 'Color',[0.85 0.2 0.2],'LineWidth',2);
        title(sprintf('Negative Linear Overlap (N=%d)\nr=%.3f, p=%.2e',numel(t_lin_vec_gm), r_lin_t, p_lin_t),...
            'FontSize',20,'FontWeight','bold');
    else
        title(sprintf('Negative Linear Overlap (N=%d)',numel(t_lin_vec_gm)),'FontSize',20,'FontWeight','bold');
    end
    xlabel('GM-GM t(Age)','FontSize',18,'FontWeight','bold');
    ylabel('GWM-HFN t(Age)','FontSize',18,'FontWeight','bold');
    set(gca,'FontSize',15,'LineWidth',1.1,'Box','on','TickDir','out','FontWeight','bold');
    axis square;
end

% Panel B: Inverted-U quadratic
nexttile;
if isempty(t_quad_vec_gm)
    text(0.5,0.5,'No inverted-U quadratic overlap edges','HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
    axis off
else
    scatter(t_quad_vec_gm, t_quad_vec_gwm, 40, 'MarkerFaceColor',[0.10 0.55 0.40], ...
        'MarkerEdgeColor','k','LineWidth',0.4,'MarkerFaceAlpha',0.75); hold on;
    mn2 = min([t_quad_vec_gm; t_quad_vec_gwm]); mx2 = max([t_quad_vec_gm; t_quad_vec_gwm]);
    if mx2==mn2, mn2 = mn2 - 0.5; mx2 = mx2 + 0.5; end
    pad2 = 0.05*(mx2-mn2);
    xlim([mn2-pad2 mx2+pad2]); ylim([mn2-pad2 mx2+pad2]);
    plot([mn2 mx2],[mn2 mx2],'k--','LineWidth',1.1);
    if numel(t_quad_vec_gm) >= 2
        mdlQ = fitlm(t_quad_vec_gm, t_quad_vec_gwm);
        xf2 = linspace(xlim(gca)*[1;0], xlim(gca)*[0;1], 200);
        y_fit2 = predict(mdlQ, xf2');
        plot(xf2, y_fit2,'Color',[0.85 0.2 0.2],'LineWidth',2);
        title(sprintf('Inverted-U Quadratic Overlap (N=%d)\nr=%.3f, p=%.2e', ...
            numel(t_quad_vec_gm), r_quad_t, p_quad_t),...
            'FontSize',20,'FontWeight','bold');
    else
        title(sprintf('Inverted-U Quadratic Overlap (N=%d)', numel(t_quad_vec_gm)),...
            'FontSize',20,'FontWeight','bold');
    end
    xlabel('GM-GM t(Age^2)','FontSize',18,'FontWeight','bold');
    ylabel('GWM-HFN t(Age^2)','FontSize',18,'FontWeight','bold');
    set(gca,'FontSize',15,'LineWidth',1.1,'Box','on','TickDir','out','FontWeight','bold');
    axis square;
end

if save_fig
    print(fullfile(fig_dir,'Figure_Overlap_tScatter_NegLinear_InvertedU_Combined.tiff'),'-dtiff','-r300');
    print(fullfile(fig_dir,'Figure_Overlap_tScatter_NegLinear_InvertedU_Combined.svg'),'-dsvg');
end

tStats = struct();
tStats.alpha_fdr                = alpha_fdr;
tStats.Linear.N_negativeOverlap = numel(t_lin_vec_gm);
tStats.Linear.r_raw             = r_lin_t;
tStats.Linear.p_raw             = p_lin_t;
tStats.Linear.r_z               = r_lin_z;
tStats.Linear.p_z               = p_lin_z;
tStats.Quadratic.N_invertedU    = numel(t_quad_vec_gm);
tStats.Quadratic.r_raw          = r_quad_t;
tStats.Quadratic.p_raw          = p_quad_t;
tStats.Quadratic.r_z            = r_quad_z;
tStats.Quadratic.p_z            = p_quad_z;
tStats.useSizeEncoding          = useSizeEncoding;

save(fullfile(fig_dir,'Overlap_tStats_NegLinear_InvertedU_NoColorbar.mat'), ...
    'tStats','idx_i_lin','idx_j_lin','idx_i_q','idx_j_q');

fprintf('\nCompleted (negative linear & inverted-U quadratic, no colorbar).\n');
fprintf('Negative linear overlap: N=%d, r=%.3f (p=%.2e)\n', ...
    tStats.Linear.N_negativeOverlap, tStats.Linear.r_raw, tStats.Linear.p_raw);
fprintf('Inverted-U quadratic overlap: N=%d, r=%.3f (p=%.2e)\n', ...
    tStats.Quadratic.N_invertedU, tStats.Quadratic.r_raw, tStats.Quadratic.p_raw);



%% ===================== (10) Network-Level Comparison =====================
net_index_raw = readmatrix(network_idx_file);
net_index_raw = net_index_raw(:)';
assert(numel(net_index_raw)==90,'Network index length must be 90.');
[~,~,net_idx_1K] = unique(net_index_raw);
net_index = net_idx_1K;
K = numel(unique(net_index));

Den = zeros(K,K);
for a = 1:K
    nodes_a = find(net_index==a);
    na = numel(nodes_a);
    Den(a,a) = na*(na-1)/2;
    for b = a+1:K
        nodes_b = find(net_index==b);
        nb = numel(nodes_b);
        Den(a,b) = na*nb;
        Den(b,a) = Den(a,b);
    end
end

Count_GWM     = zeros(K,K);
Count_GM      = zeros(K,K);
Count_Overlap = zeros(K,K);

for kEdge = 1:n_edges
    i = idx_i(kEdge); j = idx_j(kEdge);
    a = net_index(i); b = net_index(j);
    if sig_gwm_LT(i,j)
        Count_GWM(a,b) = Count_GWM(a,b)+1;
        if a~=b, Count_GWM(b,a)=Count_GWM(b,a)+1; end
    end
    if sig_gm_LT(i,j)
        Count_GM(a,b) = Count_GM(a,b)+1;
        if a~=b, Count_GM(b,a)=Count_GM(b,a)+1; end
    end
    if overlap_sig_LT(i,j)
        Count_Overlap(a,b) = Count_Overlap(a,b)+1;
        if a~=b, Count_Overlap(b,a)=Count_Overlap(b,a)+1; end
    end
end

Prop_GWM     = Count_GWM ./ max(Den,1);
Prop_GM      = Count_GM  ./ max(Den,1);
Prop_Overlap = Count_Overlap ./ max(Den,1);

NetA=[]; NetB=[]; DenVec=[]; CountGWM=[]; CountGM=[]; CountBoth=[]; PropGWM=[]; PropGM_vec=[]; PropBoth=[]; PairLabel={};
for a=1:K
    for b=a:K
        NetA(end+1,1)=a; %#ok<AGROW>
        NetB(end+1,1)=b; %#ok<AGROW>
        DenVec(end+1,1)=Den(a,b);
        CountGWM(end+1,1)=Count_GWM(a,b);
        CountGM(end+1,1)=Count_GM(a,b);
        CountBoth(end+1,1)=Count_Overlap(a,b);
        PropGWM(end+1,1)=Prop_GWM(a,b);
        PropGM_vec(end+1,1)=Prop_GM(a,b);
        PropBoth(end+1,1)=Prop_Overlap(a,b);
        PairLabel{end+1,1}=sprintf('%d-%d',a,b); %#ok<AGROW>
    end
end

NetworkPairTable = table(NetA, NetB, DenVec, CountGWM, CountGM, CountBoth, ...
    PropGWM, PropGM_vec, PropBoth, PairLabel, ...
    'VariableNames', {'NetA','NetB','Denominator','Count_GWMHFN','Count_GMGM','Count_Overlap',...
    'Prop_GWMHFN','Prop_GMGM','Prop_Overlap','PairLabel'});
writetable(NetworkPairTable, fullfile(out_dir,'NetworkPair_SignificanceComparison.csv'));

figure('Color','w','Position',[100 100 900 300]);
subplot(1,3,1); imagesc(Prop_GWM); axis square; colorbar; title('Prop GWM-HFN');
subplot(1,3,2); imagesc(Prop_GM); axis square; colorbar; title('Prop GM-GM');
subplot(1,3,3); imagesc(Prop_GWM - Prop_GM); axis square; colorbar; title('Difference (GWM-HFN - GM-GM)');
colormap(parula);
% if do_save_figures
%     print(fullfile(fig_dir,'NetworkPair_PropComparison.tiff'),'-dtiff','-r300');
% end

%% ===================== (11) Subject-level Method Agreement vs Age =====================
raw_GWM_file = fullfile(data_root,'FCMData_GWMHFN.mat');
raw_GM_file  = fullfile(data_root,'FCMData_GMGM.mat');
demo_file    = fullfile(data_root,'Anlysis_DemographicInfo.csv');

if exist(raw_GWM_file,'file') && exist(raw_GM_file,'file')
    SW_GWM = load(raw_GWM_file);
    SW_GM  = load(raw_GM_file);
    if isfield(SW_GWM,'fcm_GWG_cstd_prod') && isfield(SW_GM,'Z_GM')
        Z_GWMHFN_raw = SW_GWM.fcm_GWG_cstd_prod;
        Z_GM_raw     = SW_GM.Z_GM;
        if isequal(size(Z_GWMHFN_raw), size(Z_GM_raw))
            Nsubj = size(Z_GM_raw,3);
            Tdemo = readtable(demo_file,'TextType','string');
            Age_vec = double(Tdemo.Age);
            if numel(Age_vec)~=Nsubj
                warning('Demographics length mismatch.');
                Age_vec = Age_vec(1:min(end,Nsubj));
            end
            gm_mat  = zeros(Nsubj, n_edges);
            gwm_mat = zeros(Nsubj, n_edges);
            for e = 1:n_edges
                i = idx_i(e); j = idx_j(e);
                gm_mat(:,e)  = squeeze(Z_GM_raw(i,j,:));
                gwm_mat(:,e) = squeeze(Z_GWMHFN_raw(i,j,:));
            end
            subject_correlations = zeros(Nsubj,1);
            subject_p = zeros(Nsubj,1);
            for s=1:Nsubj
                [r_s,p_s] = corr(gm_mat(s,:)', gwm_mat(s,:)','type','Pearson','rows','pairwise');
                subject_correlations(s)=r_s;
                subject_p(s)=p_s;
            end
            [r_age,p_age] = corr(Age_vec, subject_correlations,'type','Pearson','rows','pairwise');
            Age_z = (Age_vec - mean(Age_vec,'omitnan'))/std(Age_vec,'omitnan');
            mdl_subj = fitlm(Age_z, subject_correlations);
            beta_age_subj   = mdl_subj.Coefficients.Estimate(2);
            p_beta_age_subj = mdl_subj.Coefficients.pValue(2);

            figure('Position', [100, 100, 800, 600],'Color','w');
            pos_scatter = [0.15 0.15 0.6 0.7];
            pos_right   = [0.75 0.15 0.15 0.7];
            main_ax = axes('Position', pos_scatter);
            scatter(Age_vec, subject_correlations, 55,'filled','MarkerFaceColor',[0.2 0.4 0.8], ...
                'MarkerFaceAlpha',0.65,'MarkerEdgeColor','none'); hold on;
            x_range = linspace(18,80,100);
            y_fit   = mdl_subj.Coefficients.Estimate(1) + mdl_subj.Coefficients.Estimate(2)*((x_range-mean(Age_vec))/std(Age_vec));
            plot(x_range, y_fit,'LineWidth',2.5,'Color',[0.8 0.2 0.2]);
            xlabel('Age (years)','FontSize',20,'FontWeight','bold');
            ylabel('Between-Method Correlation','FontSize',16,'FontWeight','bold');
            set(main_ax,'FontSize',16,'FontWeight','bold','LineWidth',1.5,'Box','on',...
                'TickDir','out','TickLength',[0.015 0.015]);
            grid on;
            y_lim = ylim;

            right_ax = axes('Position', pos_right);
            [counts, edges_hist] = histcounts(subject_correlations, 20,'Normalization','probability');
            centers = (edges_hist(1:end-1)+edges_hist(2:end))/2;
            barh(centers, counts, 1,'FaceColor',[0.2 0.4 0.8],'EdgeColor','none','FaceAlpha',0.35); hold on;
            [f_d, xi_d] = ksdensity(subject_correlations);
            f_d = f_d * (max(counts)/max(f_d));
            plot(f_d, xi_d,'LineWidth',2,'Color',[0.1 0.2 0.5]);
            set(right_ax,'YLim',y_lim,'XTickLabel',[],'YTickLabel',[],'Color','none','Box','off',...
                'TickDir','out','LineWidth',1.0);

            axes(main_ax);
            hTitle = title('Age-Related Changes in Method Agreement','FontSize',22,'FontWeight','bold');
            set(hTitle,'Units','normalized','Position',[0.5 1.07 0]);
            stats_str = sprintf('r = %.3f, p = %.1e', r_age, p_age);
            text(0.05,0.95,stats_str,'Units','normalized','FontSize',18,'FontWeight','bold',...
                'BackgroundColor',[1 1 1 0.8],'EdgeColor',[0.5 0.5 0.5],'Margin',5);

            if do_save_figures
                print(fullfile(fig_dir,'SubjectAgreement_Age.tiff'),'-dtiff','-r300');
                saveas(gcf, fullfile(fig_dir,'SubjectAgreement_Age.fig'));
                print(fullfile(fig_dir,'SubjectAgreement_Age.svg'),'-dsvg');
            end

            SubjectAgreement = table((1:Nsubj)', Age_vec(:), subject_correlations(:), subject_p(:), ...
                'VariableNames', {'SubjectIdx','Age','MethodCorrelation','P_value'});
            SubjectAgreementStats = table(beta_age_subj, p_beta_age_subj, r_age, p_age, ...
                'VariableNames', {'Beta_Age_z','Beta_Age_p','r_age','r_age_p'});
            writetable(SubjectAgreement,      fullfile(out_dir,'SubjectLevel_MethodAgreement.csv'));
            writetable(SubjectAgreementStats, fullfile(out_dir,'SubjectLevel_MethodAgreement_Summary.csv'));
        else
            warning('Raw matrix dimension mismatch. Skipping subject-level agreement.');
        end
    else
        warning('Expected raw matrices not found in loaded files.');
    end
else
    warning('Raw connectivity files not found. Skipping subject-level agreement.');
end

%% ===================== (12) Consolidate & Save Master MAT =====================
Comparison = struct();
Comparison.Params              = struct('alpha_fdr',alpha_fdr);
Comparison.EdgeOverlapSummary = EdgeOverlapSummary;
Comparison.DetailedOverlap    = DetailedOverlap;
Comparison.NetworkPairTable   = NetworkPairTable;
Comparison.CountMatrices      = struct('Count_GWMHFN',Count_GWM,'Count_GMGM',Count_GM,...
    'Count_Overlap',Count_Overlap,'Denominator',Den);
if exist('SubjectAgreement','var')
    Comparison.SubjectAgreement      = SubjectAgreement;
    Comparison.SubjectAgreementStats = SubjectAgreementStats;
end
Comparison.tStats = tStats; %#ok<NODEF>

save(fullfile(out_dir,'Comparison_Age_GWMHFN_vs_GMGM.mat'),'Comparison','-v7.3');
fprintf('\nAll comparison outputs saved in: %s\n', out_dir);