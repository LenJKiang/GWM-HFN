% Analysis_CorrelationWithAge_GMGM_latest.m
% Age effects for GM-GM connectivity in SALD with formula-free fitlm calls
% - Robust data import (Sex/Hand safe)
% - Global mean: linear vs quadratic (via explicit X matrices, no formula strings)
% - Edgewise model selection by AIC (linear vs quadratic), BH-FDR
% - Network-level proportions (within vs between) with two-proportion logic
% - Circos link files (linear and quadratic)
% - Structured outputs for downstream comparison
%
% Author: Zeqiang Linli
% Date: 2025-08-13

clear; clc; close all;

%% ===================== User configuration =====================
data_root = 'E:\Neuroimage\MyProject\GMWM_Network\Data\SALD';
mat_file  = fullfile(data_root, 'FCMData_GMGM.mat');             % contains Z_GM (90x90xN), global_GM (Nx1)
demo_csv  = fullfile(data_root, 'Anlysis_DemographicInfo.csv');  % demographics

% Network index (AAL90 7-network assignment)
network_index_csv = 'E:\Neuroimage\MyProject\GMWM_Network\Data\AAL_7networkIndex.csv';
network_names = {'VN','SMN','AN','LB','FPN','DMN','BGN'};

% Circos genes info (AAL)
genes_info_file = 'E:\Neuroimage\Software\circos-0.69-6\AAL\genes.labelsAAL.links.txt';

% Output directory
out_dir   = 'E:\Neuroimage\MyProject\GMWM_Network\Code_R1\step4_association_age';
fig_dir   = fullfile('E:\Neuroimage\MyProject\GMWM_Network\Code_R1\visualization\Figure6_Age');
circos_dir= fullfile(out_dir, 'circos');
if ~exist(out_dir,'dir');    mkdir(out_dir);    end
if ~exist(fig_dir,'dir');    mkdir(fig_dir);    end
if ~exist(circos_dir,'dir'); mkdir(circos_dir); end

% FDR alpha for edgewise significance and visualization
alpha_fdr = 0.05;

%% ===================== Load connectivity and demographics =====================
% Load matrices
S = load(mat_file);
if ~isfield(S, 'Z_GM') || ~isfield(S, 'global_GM')
    error('Expected variables Z_GM (90x90xN) and global_GM (Nx1) in %s', mat_file);
end
Z_GM      = S.Z_GM;             % 90x90xN
Global_GM = S.global_GM(:);     % Nx1

[n1, n2, N] = size(Z_GM);
if n1~=90 || n2~=90, error('Z_GM must be 90x90xN'); end
if numel(Global_GM)~=N, error('Global vector length does not match N'); end

% Load demographics safely (Sex/Hand can be string/categorical/cell)
D = readtable(demo_csv, 'TextType','string');
if height(D)~=N
    warning('Demographic rows (%d) != subjects (%d). Ensure ordering matches subjects in matrices.', height(D), N);
end

% Build covariates
Age   = double(D.Age);
FD    = double(D.meanFD);
SexM  = double(ismember(upper(string(D.Sex)), ["M","MALE"]));
HandR = double(ismember(upper(string(D.Handedness)), ["RIGHT","R"]));
% Optional: mark missing as NaN (fitlm will omit rows)
SexM(ismissing(D.Sex)) = NaN;
HandR(ismissing(D.Handedness)) = NaN;

% Z-score continuous covariates
AgeZ  = zscore(Age);
FDZ   = zscore(FD);
AgeZ2 = AgeZ.^2; % quadratic term

%% ===================== Global mean analysis (no formula strings) =====================
% Build design matrices
Xg_lin  = [AgeZ,     SexM, HandR, FDZ];
Xg_quad = [AgeZ,AgeZ2,SexM, HandR, FDZ];
Yg = Global_GM;

% Mask valid rows (drop any NaNs)
mask_lin  = all(~isnan([Xg_lin,  Yg]),  2);
mask_quad = all(~isnan([Xg_quad, Yg]),  2);

% Provide VarNames so coefficients have readable labels
varnames_lin  = {'AgeZ','SexM','HandR','FDZ','Global'};
varnames_quad = {'AgeZ','AgeZ2','SexM','HandR','FDZ','Global'};

mdl_lin  = fitlm(Xg_lin(mask_lin,:),  Yg(mask_lin),  'linear', 'VarNames', varnames_lin);
mdl_quad = fitlm(Xg_quad(mask_quad,:), Yg(mask_quad), 'linear', 'VarNames', varnames_quad);

AIC_lin  = mdl_lin.ModelCriterion.AIC;
AIC_quad = mdl_quad.ModelCriterion.AIC;
if AIC_quad < AIC_lin
    best_model = 'Quadratic';
else
    best_model = 'Linear';
end

fprintf('\n=== Global mean GM-GM analysis ===\n');
fprintf('AIC (Linear)   = %.3f\n', AIC_lin);
fprintf('AIC (Quadratic)= %.3f\n', AIC_quad);
fprintf('Best model     = %s\n', best_model);
% AIC (Linear)   = -511.755
% AIC (Quadratic)= -514.532
% Best model     = Quadratic
%% ===================== Global mean scatter plot (Linear model preferred) =====================
% Draw scatter (Age vs Global GM-GM) with density, regression line, and partial correlation.
% Note:
% - Requires "scatplot" function available on MATLAB path.
% - Saves figure to fig_dir.

if strcmpi(best_model, 'Quadratic')
    Y = Global_GM; 

    valid_mask = ~isnan(Age) & ~isnan(Y) & ~isnan(SexM) & ~isnan(HandR) & ~isnan(FD);
    Age_v  = Age(valid_mask);
    Y_v    = Y(valid_mask);
    Sex_v  = SexM(valid_mask);
    Hand_v = HandR(valid_mask);
    FD_v   = FD(valid_mask);

    Age_mean = mean(Age_v);
    Age_c  = Age_v - Age_mean;
    Age_c2 = Age_c.^2;

    tbl = table(Y_v, Age_c, Age_c2, Sex_v, Hand_v, FD_v, ...
        'VariableNames', {'Y','Age_c','Age_c2','Sex','Hand','FD'});
    mdl_q = fitlm(tbl, 'Y ~ Age_c + Age_c2 + Sex + Hand + FD');

    % mdl_lin = fitlm(tbl, 'Y ~ Age_c + Sex + Hand + FD');

    figure('Position', [100, 100, 800, 600], 'Color','w');

    % Density scatter
    out = scatplot(Age_v, Y_v, [], [], [], [], [], 25);
    colorbar('off'); hold on;

    x_plot = linspace(min(Age_v), max(Age_v), 200)';
    Age_c_plot  = x_plot - Age_mean;
    Age_c2_plot = Age_c_plot.^2;

    sex_mode  = mode(Sex_v);   % 若 Sex_v 是 0/1
    hand_mode = mode(Hand_v);
    fd_mean   = mean(FD_v);

    newTbl = table(Age_c_plot, Age_c2_plot, ...
        repmat(sex_mode, size(x_plot)), ...
        repmat(hand_mode, size(x_plot)), ...
        repmat(fd_mean, size(x_plot)), ...
        'VariableNames', {'Age_c','Age_c2','Sex','Hand','FD'});

    [y_fit, y_ci] = predict(mdl_q, newTbl);  % y_ci: n x 2 (lower upper)

    fill([x_plot; flipud(x_plot)], [y_ci(:,1); flipud(y_ci(:,2))], ...
        [0.8 0.2 0.2], 'FaceAlpha', 0.15, 'EdgeColor','none');

    plot(x_plot, y_fit, 'LineWidth', 2.8, 'Color', [0.8 0.2 0.2]);

    b1 = mdl_q.Coefficients.Estimate(strcmp(mdl_q.Coefficients.Row,'Age_c'));
    b2 = mdl_q.Coefficients.Estimate(strcmp(mdl_q.Coefficients.Row,'Age_c2'));
    if b2 ~= 0
        Age_vertex = Age_mean - b1 / (2*b2);
        if Age_vertex >= min(Age_v) && Age_vertex <= max(Age_v)
            y_vertex = predict(mdl_q, table(Age_vertex - Age_mean, (Age_vertex - Age_mean)^2, sex_mode, hand_mode, fd_mean, ...
                'VariableNames', {'Age_c','Age_c2','Sex','Hand','FD'}));
            plot(Age_vertex, y_vertex, 'ko', 'MarkerFaceColor','k', 'MarkerSize',6);
            text(Age_vertex, y_vertex, sprintf('  Peak %.1f y', Age_vertex), ...
                'FontSize', 16, 'FontWeight','bold', 'Color','k', ...
                'VerticalAlignment','bottom');
        end
    end

    xlim([min(Age_v)-1, max(Age_v)+1]);
    y1 = prctile(Y_v,1); y2 = prctile(Y_v,99);
    pad = 0.05 * max(eps, y2 - y1);
    ylim([y1 - pad, y2 + pad]);

    t_quad = mdl_q.Coefficients.tStat(strcmp(mdl_q.Coefficients.Row,'Age_c2'));
    p_quad = mdl_q.Coefficients.pValue(strcmp(mdl_q.Coefficients.Row,'Age_c2'));
    df     = mdl_q.DFE;
    partial_R2_quad = t_quad^2 / (t_quad^2 + df);



    if p_quad < 1e-3
        p_txt = 'p < 1e-3';
    else
        p_txt = sprintf('p = %.3f', p_quad);
    end
    annot_str = sprintf('Quadratic term (Age^2): %s', p_txt);
    text(0.05, 0.92, annot_str, 'Units','normalized', ...
        'FontSize', 20, 'FontWeight','bold', ...
        'BackgroundColor',[1 1 1 0.80], ...
        'EdgeColor',[0.5 0.5 0.5], 'Margin',6);
    set(gca,'FontSize',16,'FontWeight','bold','LineWidth',1.5,...
        'Box','on','TickDir','out','TickLength',[0.01 0.01],...
        'GridColor',[0.8 0.8 0.8],'GridAlpha',0.3);
    title('Correlation of Mean GM-GM and Age', ...
        'FontSize', 25, 'FontWeight','bold', ...
        'Units','normalized','Position',[0.5,1.02]);
    xlabel('Age (years)', 'FontSize', 23, 'FontWeight','bold');
    ylabel('Mean GM-GM',   'FontSize', 23, 'FontWeight','bold');


    grid on;

    set(gcf,'Position',[100,100,700,600]);

    out_name = fullfile(fig_dir, 'Scatterplot_meanGMGM_Age_Quadratic');
    print([out_name '.tiff'],'-dtiff','-r300');
    saveas(gcf, [out_name '.fig']);
    print([out_name '.svg'],'-dsvg');

    clear Age_v Y_v Sex_v Hand_v FD_v Age_c Age_c2 tbl mdl_q ...
          x_plot Age_c_plot Age_c2_plot newTbl y_fit y_ci ...
          b1 b2 Age_vertex y_vertex t_quad p_quad partial_R2_quad annot_str p_txt pad
end

%% ===================== Edgewise analysis (AIC: linear vs quadratic; no formula strings) =====================
% Prepare lower triangle indices
maskLT = tril(true(90), -1);
[idx_i, idx_j] = find(maskLT);
n_edges = numel(idx_i);

% Pre-allocate
beta_age   = nan(90,90);   t_age   = nan(90,90);   p_age   = nan(90,90);
beta_age2  = nan(90,90);   t_age2  = nan(90,90);   p_age2  = nan(90,90);
model_used = zeros(90,90); % 1=linear, 2=quadratic
AIC_L_mat  = nan(90,90);   AIC_Q_mat= nan(90,90);
peak_age_yrs = nan(90,90); % vertex for quadratic in real years

% Reuse covariates matrices
Xe_lin  = [AgeZ,     SexM, HandR, FDZ];
Xe_quad = [AgeZ,AgeZ2,SexM, HandR, FDZ];

fprintf('\nFitting %d GM-GM edges (linear vs quadratic)...\n', n_edges);
for k = 1:n_edges
    i = idx_i(k); j = idx_j(k);
    Y = squeeze(Z_GM(i,j,:));  % Nx1

    % Linear
    maskL = all(~isnan([Xe_lin,  Y]), 2);
    mdl_L = fitlm(Xe_lin(maskL,:),  Y(maskL),  'linear', 'VarNames', {'AgeZ','SexM','HandR','FDZ','Y'});

    % Quadratic
    maskQ = all(~isnan([Xe_quad, Y]), 2);
    mdl_Q = fitlm(Xe_quad(maskQ,:), Y(maskQ), 'linear', 'VarNames', {'AgeZ','AgeZ2','SexM','HandR','FDZ','Y'});

    AIC_L = mdl_L.ModelCriterion.AIC;
    AIC_Q = mdl_Q.ModelCriterion.AIC;
    AIC_L_mat(i,j) = AIC_L;
    AIC_Q_mat(i,j) = AIC_Q;

    if AIC_Q < AIC_L
        model_used(i,j) = 2;
        % AgeZ and AgeZ^2 terms
        c1 = mdl_Q.Coefficients.Estimate('AgeZ');
        t1 = mdl_Q.Coefficients.tStat('AgeZ');
        p1 = mdl_Q.Coefficients.pValue('AgeZ');

        c2 = mdl_Q.Coefficients.Estimate('AgeZ2');
        t2 = mdl_Q.Coefficients.tStat('AgeZ2');
        p2 = mdl_Q.Coefficients.pValue('AgeZ2');

        beta_age(i,j)  = c1; t_age(i,j)  = t1; p_age(i,j)  = p1;
        beta_age2(i,j) = c2; t_age2(i,j) = t2; p_age2(i,j) = p2;

        if ~isnan(c2) && c2~=0
            peak_age_yrs(i,j) = -c1/(2*c2) * std(Age) + mean(Age);
        end
    else
        model_used(i,j) = 1;
        c1 = mdl_L.Coefficients.Estimate('AgeZ');
        t1 = mdl_L.Coefficients.tStat('AgeZ');
        p1 = mdl_L.Coefficients.pValue('AgeZ');

        beta_age(i,j) = c1; t_age(i,j) = t1; p_age(i,j) = p1;
    end
end

% Combined p for FDR
combined_p = nan(90,90);
lin_sel  = (model_used==1);
quad_sel = (model_used==2);
combined_p(lin_sel)  = p_age(lin_sel);
combined_p(quad_sel) = p_age2(quad_sel);

% Benjamini-Hochberg FDR across lower triangle
p_vec = combined_p(maskLT);
valid = ~isnan(p_vec);
m = sum(valid);
q_vec = nan(size(p_vec));
if m > 0
    pv  = p_vec(valid);
    [pv_sorted, order] = sort(pv);
    ranks = (1:m)';
    q_sorted = (m ./ ranks) .* pv_sorted;
    q_sorted = flipud(cummin(flipud(q_sorted)));
    q_sorted = min(q_sorted, 1);
    invOrder = zeros(m,1); invOrder(order) = 1:m;
    q_vec(valid) = q_sorted(invOrder);
end
q_mat = nan(90,90);
q_mat(maskLT) = q_vec;
p_mat = nan(90,90);
p_mat(maskLT) = p_vec;

% Significant masks (q < alpha)
sig_mask_all = false(90,90);
sig_mask_all(maskLT) = q_vec < alpha_fdr;
sig_mask_all = sig_mask_all | sig_mask_all.'; % symmetric for convenience

sig_mask_linear    = (model_used==1) & (q_mat<alpha_fdr);
sig_mask_quadratic = (model_used==2) & (q_mat<alpha_fdr);

fprintf('Edgewise results (GM-GM): %d/%d significant edges at q<%.2f.\n', ...
    sum((q_vec<alpha_fdr), 'omitnan'), nnz(maskLT), alpha_fdr);
% Fitting 4005 GM-GM edges (linear vs quadratic)...
% Edgewise results (GM-GM): 1550/4005 significant edges at q<0.05.
%% ===================== Summary of significant edge pattern types (Linear / Quadratic) =====================

% Lower-triangle significant edge mask (avoid double counting)
sig_edges_LT = (q_mat < alpha_fdr) & maskLT;

% Masks separating linear- and quadratic-selected significant edges
sig_lin_edges_LT  = sig_edges_LT & (model_used == 1);
sig_quad_edges_LT = sig_edges_LT & (model_used == 2);

% Counts: linear positive / linear negative
n_linear_pos = sum(beta_age(sig_lin_edges_LT)  > 0, 'omitnan');
n_linear_neg = sum(beta_age(sig_lin_edges_LT)  < 0, 'omitnan');

% Counts: quadratic U-shape (Age^2 coefficient > 0) / inverted U (Age^2 coefficient < 0)
n_quad_u     = sum(beta_age2(sig_quad_edges_LT) > 0, 'omitnan');  % convex
n_quad_inv_u = sum(beta_age2(sig_quad_edges_LT) < 0, 'omitnan');  % concave

% Totals
n_sig_total     = sum(sig_edges_LT(:));
n_sig_linear    = sum(sig_lin_edges_LT(:));
n_sig_quadratic = sum(sig_quad_edges_LT(:));

% Optional: restrict to quadratic edges whose peak (vertex) is within the empirical age range
if exist('Age','var') && ~isempty(Age)
    age_min = min(Age);
    age_max = max(Age);
    quad_sig_mask_full = (model_used == 2) & sig_edges_LT;
    peak_in_range_mask = quad_sig_mask_full & ...
        (peak_age_yrs >= age_min) & (peak_age_yrs <= age_max);

    n_quad_u_inRange     = sum((beta_age2 > 0) & peak_in_range_mask, 'omitnan');
    n_quad_inv_u_inRange = sum((beta_age2 < 0) & peak_in_range_mask, 'omitnan');
else
    age_min = NaN; age_max = NaN;
    n_quad_u_inRange = NaN;
    n_quad_inv_u_inRange = NaN;
end

% Percentages relative to all significant edges
pct_total = @(x) 100 * x / max(1, n_sig_total);
pct_linear_pos_total = pct_total(n_linear_pos);
pct_linear_neg_total = pct_total(n_linear_neg);
pct_quad_u_total     = pct_total(n_quad_u);
pct_quad_inv_u_total = pct_total(n_quad_inv_u);

fprintf('\n=== GM-GM Significant Edge Pattern Summary (q < %.2f) ===\n', alpha_fdr);
fprintf('Total significant edges: %d\n', n_sig_total);
fprintf('  Linear-selected:   %d (%.1f%% of significant)\n', ...
    n_sig_linear, 100*n_sig_linear/max(1,n_sig_total));
fprintf('    Positive linear: %d (%.1f%% of significant, %.1f%% of linear)\n', ...
    n_linear_pos, pct_linear_pos_total, 100*n_linear_pos/max(1,n_sig_linear));
fprintf('    Negative linear: %d (%.1f%% of significant, %.1f%% of linear)\n', ...
    n_linear_neg, pct_linear_neg_total, 100*n_linear_neg/max(1,n_sig_linear));
fprintf('  Quadratic-selected:%d (%.1f%% of significant)\n', ...
    n_sig_quadratic, 100*n_sig_quadratic/max(1,n_sig_total));
fprintf('    U-shape (Age^2 > 0):      %d (%.1f%% of significant, %.1f%% of quadratic)\n', ...
    n_quad_u, pct_quad_u_total, 100*n_quad_u/max(1,n_sig_quadratic));
fprintf('    Inverted U (Age^2 < 0):   %d (%.1f%% of significant, %.1f%% of quadratic)\n', ...
    n_quad_inv_u, pct_quad_inv_u_total, 100*n_quad_inv_u/max(1,n_sig_quadratic));

if ~isnan(age_min)
    fprintf('    (Vertex within observed age range [%g, %g]):\n', age_min, age_max);
    fprintf('       U-shape in-range:      %d\n', n_quad_u_inRange);
    fprintf('       Inverted U in-range:   %d\n', n_quad_inv_u_inRange);
end

% Total significant edges: 1550
%   Linear-selected:   656 (42.3% of significant)
%     Positive linear: 157 (10.1% of significant, 23.9% of linear)
%     Negative linear: 499 (32.2% of significant, 76.1% of linear)
%   Quadratic-selected:894 (57.7% of significant)
%     U-shape (Age^2 > 0):      5 (0.3% of significant, 0.6% of quadratic)
%     Inverted U (Age^2 < 0):   889 (57.4% of significant, 99.4% of quadratic)

% Append to AgeEffects structure if it already exists
if exist('AgeEffects','var') && isfield(AgeEffects,'Edgewise')
    AgeEffects.Edgewise.PatternSummary = struct( ...
        'n_sig_total', n_sig_total, ...
        'n_sig_linear', n_sig_linear, ...
        'n_linear_pos', n_linear_pos, ...
        'n_linear_neg', n_linear_neg, ...
        'n_sig_quadratic', n_sig_quadratic, ...
        'n_quad_u', n_quad_u, ...
        'n_quad_inv_u', n_quad_inv_u, ...
        'n_quad_u_inRange', n_quad_u_inRange, ...
        'n_quad_inv_u_inRange', n_quad_inv_u_inRange, ...
        'pct_linear_pos_total', pct_linear_pos_total, ...
        'pct_linear_neg_total', pct_linear_neg_total, ...
        'pct_quad_u_total', pct_quad_u_total, ...
        'pct_quad_inv_u_total', pct_quad_inv_u_total);
end

%% ===================== Save edgewise table (CSV) =====================
% Load network index
net_index = readmatrix(network_index_csv);
net_index = net_index(:)';   % force row
if numel(net_index)~=90, error('Network index must have 90 elements'); end

edge_rows = idx_i; edge_cols = idx_j;
T_edge = table(edge_rows, edge_cols, ...
    net_index(edge_rows)', net_index(edge_cols)', ...
    model_used(maskLT), ...
    beta_age(maskLT),   t_age(maskLT),   p_age(maskLT), ...
    beta_age2(maskLT),  t_age2(maskLT),  p_age2(maskLT), ...
    q_mat(maskLT),      peak_age_yrs(maskLT), ...
    'VariableNames', {'Node1','Node2','Net1','Net2','ModelUsed', ...
                      'Beta_Age','t_Age','p_Age', ...
                      'Beta_Age2','t_Age2','p_Age2', ...
                      'q_value','PeakAge_years'});

csv_edge = fullfile(out_dir, 'Edgewise_AgeEffects_GMGM.csv');
writetable(T_edge, csv_edge);
fprintf('Saved edgewise results to: %s\n', csv_edge);

%% ===================== Build significant matrices (for visualization and circos) =====================
matrix_age  = zeros(90,90);
matrix_age2 = zeros(90,90);

sig_lin_idx  = find((model_used==1) & (q_mat<alpha_fdr) & maskLT);
sig_quad_idx = find((model_used==2) & (q_mat<alpha_fdr) & maskLT);
[rL,cL] = ind2sub([90,90], sig_lin_idx);
[rQ,cQ] = ind2sub([90,90], sig_quad_idx);

for u = 1:numel(sig_lin_idx)
    matrix_age(rL(u), cL(u)) = beta_age(rL(u), cL(u));
end
for u = 1:numel(sig_quad_idx)
    matrix_age2(rQ(u), cQ(u)) = beta_age2(rQ(u), cQ(u));
end

save(fullfile(out_dir,'Matrix_Age_GMGM.mat'),  'matrix_age');
save(fullfile(out_dir,'Matrix_Age2_GMGM.mat'), 'matrix_age2');

% Quick visuals
f1 = figure('Color','w'); imagesc(matrix_age); axis square; colorbar;
title('GM-GM: Significant Age (linear) coefficients, q<0.05');
print(f1, fullfile(fig_dir,'Matrix_Age_linear_GMGM.tiff'), '-dtiff','-r300');
f2 = figure('Color','w'); imagesc(matrix_age2); axis square; colorbar;
title('GM-GM: Significant Age^2 (quadratic) coefficients, q<0.05');
print(f2, fullfile(fig_dir,'Matrix_Age_quadratic_GMGM.tiff'), '-dtiff','-r300');

%% ===================== Network-level proportions (within vs between) =====================
% Numerator_within: # significant within-network edges (count once; LT)
% Numerator_between: # significant between-network edges (count once; LT)
% Denominator_within (global): sum_k C(n_k, 2)
% Denominator_between (global): 4005 - Denominator_within (C(90,2)=4005)
% Per-network:
%   Within_den(k)  = C(n_k, 2)
%   Between_den(k) = n_k * (90 - n_k)
% Proportion = Numerator / Denominator

num_nodes = 90;
% Significant masks restricted to lower triangle (model-selected)
sig_lin_LT  = (model_used==1) & (q_mat<alpha_fdr) & maskLT;
sig_quad_LT = (model_used==2) & (q_mat<alpha_fdr) & maskLT;

% Ensure network_index is 1..K
[~, ~, net_idx_12K] = unique(net_index);
network_index = net_idx_12K;
K = numel(unique(network_index));

% Node counts per network
n_k = zeros(K,1);
for k = 1:K
    n_k(k) = sum(network_index==k);
end

% Global denominators
Den_within_global  = sum(n_k.*(n_k-1)/2);
Total_edges        = num_nodes*(num_nodes-1)/2;   % 4005
Den_between_global = Total_edges - Den_within_global;

% Global within/between masks (LT)
within_mask_LT = false(num_nodes);
for a=1:K
    nodes_a = find(network_index==a);
    tmp = false(num_nodes); tmp(nodes_a, nodes_a) = true;
    within_mask_LT = within_mask_LT | (tmp & maskLT);
end
between_mask_LT = maskLT & ~within_mask_LT;

% Global numerators
Num_within_lin   = nnz(sig_lin_LT  & within_mask_LT);
Num_between_lin  = nnz(sig_lin_LT  & between_mask_LT);
Num_within_quad  = nnz(sig_quad_LT & within_mask_LT);
Num_between_quad = nnz(sig_quad_LT & between_mask_LT);

% Global proportions
Prop_within_lin   = Num_within_lin   / max(1, Den_within_global);
Prop_between_lin  = Num_between_lin  / max(1, Den_between_global);
Prop_within_quad  = Num_within_quad  / max(1, Den_within_global);
Prop_between_quad = Num_between_quad / max(1, Den_between_global);

fprintf('\n=== GM-GM Global proportions (q < %.2f) ===\n', alpha_fdr);
fprintf('Linear:   within = %.4f, between = %.4f\n',  Prop_within_lin,  Prop_between_lin);
fprintf('Quadratic:within = %.4f, between = %.4f\n',  Prop_within_quad, Prop_between_quad);
% Linear:   within = 0.0908, between = 0.1760
% Quadratic:within = 0.4380, between = 0.1874

% Per-network numerators
linear_within_network    = zeros(K,1);
quadratic_within_network = zeros(K,1);
linear_between_pair      = zeros(K,K);
quadratic_between_pair   = zeros(K,K);

for a=1:K
    nodes_a = find(network_index==a);
    % Within
    tmp_within = false(num_nodes); tmp_within(nodes_a, nodes_a) = true;
    linear_within_network(a)    = nnz(sig_lin_LT  & tmp_within);
    quadratic_within_network(a) = nnz(sig_quad_LT & tmp_within);

    % Between to all b>a
    for b=a+1:K
        nodes_b = find(network_index==b);
        tmp_between = false(num_nodes); tmp_between(nodes_a, nodes_b) = true;
        cntL = nnz(sig_lin_LT  & tmp_between);
        cntQ = nnz(sig_quad_LT & tmp_between);
        linear_between_pair(a,b)    = cntL; linear_between_pair(b,a)    = cntL;
        quadratic_between_pair(a,b) = cntQ; quadratic_between_pair(b,a) = cntQ;
    end
end

% Per-network denominators
den_within_net  = n_k.*(n_k-1)/2;        % C(n_k,2)
den_between_net = n_k.*(num_nodes - n_k);% n_k*(90-n_k)

% Per-network proportions
prop_linear_within   = linear_within_network              ./ max(den_within_net,  1);
prop_linear_between  = sum(linear_between_pair, 2)        ./ max(den_between_net, 1);
prop_quad_within     = quadratic_within_network           ./ max(den_within_net,  1);
prop_quad_between    = sum(quadratic_between_pair, 2)     ./ max(den_between_net, 1);

% Replace any NaNs by zeros
prop_linear_within(~isfinite(prop_linear_within))   = 0;
prop_linear_between(~isfinite(prop_linear_between)) = 0;
prop_quad_within(~isfinite(prop_quad_within))       = 0;
prop_quad_between(~isfinite(prop_quad_between))     = 0;

%% --------- Figure A: Global within vs between (Linear | Quadratic) ----------
figure('Position', [80, 80, 1100, 460]); set(gcf,'Color','white');

% Linear
subplot(1,2,1);
hA = bar([Prop_within_lin, Prop_between_lin], 'grouped');
set(gca, 'FontSize', 13, 'FontName', 'Arial', 'LineWidth', 1.5, 'FontWeight', 'bold');
set(gca, 'XTickLabel', {'Within-Network','Between-Network'});
ylabel('Proportion of Significant Connections', 'FontSize', 14, 'FontWeight', 'bold');
title('GM-GM: Linear Age Effects (Global)', 'FontSize', 14, 'FontWeight', 'bold');
for i = 1:length(hA)
    xData = hA(i).XData + hA(i).XOffset;
    yData = hA(i).YData;
    for j = 1:length(xData)
        text(xData(j), yData(j), sprintf('%.1f%%', 100*yData(j)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 12, 'FontWeight', 'bold');
    end
end

% Quadratic
subplot(1,2,2);
hB = bar([Prop_within_quad, Prop_between_quad], 'grouped');
set(gca, 'FontSize', 13, 'FontName', 'Arial', 'LineWidth', 1.5, 'FontWeight', 'bold');
set(gca, 'XTickLabel', {'Within-Network','Between-Network'});
ylabel('Proportion of Significant Connections', 'FontSize', 14, 'FontWeight', 'bold');
title('GM-GM: Quadratic Age Effects (Global)', 'FontSize', 14, 'FontWeight', 'bold');
for i = 1:length(hB)
    xData = hB(i).XData + hB(i).XOffset;
    yData = hB(i).YData;
    for j = 1:length(xData)
        text(xData(j), yData(j), sprintf('%.1f%%', 100*yData(j)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 12, 'FontWeight', 'bold');
    end
end

% print(fullfile(fig_dir,'Global_WithinBetween_Proportions_GMGM.tiff'), '-dtiff', '-r300');

%% --------- Figure B: Per-network clustered bars (Linear | Quadratic) ----------
figure('Position', [100, 100, 1200, 500]); set(gcf, 'Color', 'white');

% Resolve network_names if needed
if ~exist('network_names','var') || numel(network_names)~=K
    network_names = arrayfun(@(x) sprintf('Net%d', x), 1:K, 'UniformOutput', false);
end

% Linear
subplot(1,2,1);
h1 = bar([prop_linear_within, prop_linear_between], 'grouped');
set(gca, 'FontSize', 15, 'FontName', 'Arial', 'LineWidth', 1.5, 'FontWeight', 'bold');
set(gca, 'XTickLabel', network_names);
legend('Within-Network', 'Between-Network', 'Location', 'northwest', 'FontSize', 12);
title('GM-GM: Linear Age Effects by Network', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Proportion of Significant Connections', 'FontSize', 15, 'FontWeight', 'bold');
% Labels with slight outward shifts to reduce left-right overlap
xData = h1(1).XData + h1(1).XOffset - 0.1; yData = h1(1).YData;
for j = 1:length(xData)
    text(xData(j), yData(j), sprintf('%.1f%%', 100*yData(j)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 8, 'FontWeight', 'bold');
end
xData = h1(2).XData + h1(2).XOffset + 0.1; yData = h1(2).YData;
for j = 1:length(xData)
    text(xData(j), yData(j), sprintf('%.1f%%', 100*yData(j)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 8, 'FontWeight', 'bold');
end

% Quadratic
subplot(1,2,2);
h2 = bar([prop_quad_within, prop_quad_between], 'grouped');
set(gca, 'FontSize', 15, 'FontName', 'Arial', 'LineWidth', 1.5, 'FontWeight', 'bold');
set(gca, 'XTickLabel', network_names);
legend('Within-Network', 'Between-Network', 'Location', 'northwest', 'FontSize', 12);
title('GM-GM: Quadratic Age Effects by Network', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Proportion of Significant Connections', 'FontSize', 15, 'FontWeight', 'bold');
xData = h2(1).XData + h2(1).XOffset - 0.2; yData = h2(1).YData;
for j = 1:length(xData)
    text(xData(j), yData(j), sprintf('%.1f%%', 100*yData(j)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 8, 'FontWeight', 'bold');
end
xData = h2(2).XData + h2(2).XOffset + 0.2; yData = h2(2).YData;
for j = 1:length(xData)
    text(xData(j), yData(j), sprintf('%.1f%%', 100*yData(j)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 8, 'FontWeight', 'bold');
end

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
print(fullfile(fig_dir,'Within_vs_Between_ByNetwork_Proportions_GMGM.tiff'), '-dtiff', '-r300');

%% ===================== Circos link files =====================
fprintf('\n[Circos] Generating link files (natural log p-value based)...\n');

genesInfo = readtable(genes_info_file, 'FileType','text', ...
    'Delimiter','\t','ReadVariableNames',false);
genesInfo.Properties.VariableNames = {'Network','Start','End','Region'};


% --- Configuration ---
use_qmask  = true;      % true: use q_mat < alpha_fdr as significance filter; false: raw p-value filter
alpha_raw  = 0.05;      % only used if use_qmask=false
cap_score  = 50;        % upper cap for -log(p) just for safety (won't affect bin logic)
min_p_safe = 1e-320;    % replace zeros or extremely small with this to avoid Inf

% --- Precompute lower-triangle indices ---
LT = maskLT;
[idx_i_all, idx_j_all] = find(LT);

% --- Significance masks (lower triangle only) ---
if use_qmask
    sig_lin_LT  = (model_used==1) & (q_mat < alpha_fdr) & LT;
    sig_quad_LT = (model_used==2) & (q_mat < alpha_fdr) & LT;
else
    % Raw p-value criterion: linear uses p_age, quadratic uses p_age2
    sig_lin_LT  = (model_used==1) & (p_age  < alpha_raw) & LT;
    sig_quad_LT = (model_used==2) & (p_age2 < alpha_raw) & LT;
end

n_lin  = nnz(sig_lin_LT);
n_quad = nnz(sig_quad_LT);
fprintf('[Circos] Significant edges (linear)   : %d\n', n_lin);
fprintf('[Circos] Significant edges (quadratic): %d\n', n_quad);

% Early exit if nothing
if n_lin==0 && n_quad==0
    warning('[Circos] No significant edges to write. Skipping link file creation.');
    % You may return here if desired:
    % return;
end

% --- Color tables ---
lin_pos_colors = [255 234 221; 252 174 174; 255 137 137; 255  80  80; 255  50  50];
lin_neg_colors = [223 245 255; 103 198 227;  55 140 231;  83  86 255;  50  50 231];
quad_pos_colors= [240 230 255; 229 217 242; 205 193 255; 162 148 249; 126 100 190];
quad_neg_colors= [220 255 180; 169 196 108; 128 157  60;  93 135  54;  70 100  50];

% --- Binning helper (score = -log(p), natural log) ---

% thickness map for bins 1..5
thickness_map = [1 2 3 5 9];

assign_bin = @(score) ( ...
    (score < 10) * 1 + ...
    (score >= 10 & score < 15) * 2 + ...
    (score >= 15 & score < 20) * 3 + ...
    (score >= 20 & score < 25) * 4 + ...
    (score >= 25) * 5 );

%% ---------------- Linear links ----------------
if n_lin>0
    lin_inds_vec = find(sig_lin_LT);
    file_lin = fullfile(circos_dir, 'segdup_GMGM_age_linear_gradient.link.txt');
    fid = fopen(file_lin,'w');
    if fid<0, error('Cannot open linear link file: %s', file_lin); end

    for k = 1:numel(lin_inds_vec)
        lin_linear_index = lin_inds_vec(k);    % linear index into 90x90
        [i, j] = ind2sub([90,90], lin_linear_index);

        % p for linear Age term
        p_val = p_age(i,j);
        if ~isfinite(p_val) || p_val <= 0
            p_val = min_p_safe;
        end
        score = -log(p_val);
        if score>cap_score, score = cap_score; end

        bin = assign_bin(score);
        thickness_val = thickness_map(bin);

        coeff = beta_age(i,j);  % linear coefficient
        if coeff >= 0
            color = lin_pos_colors(bin,:);
        else
            color = lin_neg_colors(bin,:);
        end

        rA = genesInfo(i,:); rB = genesInfo(j,:);
        fprintf(fid, '%s %d %d %s %d %d color=%d,%d,%d,thickness=%.2f\n', ...
            rA.Network{1}, rA.Start, rA.End, ...
            rB.Network{1}, rB.Start, rB.End, ...
            color(1), color(2), color(3), thickness_val);
    end
    fclose(fid);
else
    file_lin = '(no linear edges written)';
end

%% ---------------- Quadratic links ----------------
if n_quad>0
    quad_inds_vec = find(sig_quad_LT);
    file_quad = fullfile(circos_dir, 'segdup_GMGM_age_quadratic_gradient.link.txt');
    fid = fopen(file_quad,'w');
    if fid<0, error('Cannot open quadratic link file: %s', file_quad); end

    for k = 1:numel(quad_inds_vec)
        quad_linear_index = quad_inds_vec(k);
        [i, j] = ind2sub([90,90], quad_linear_index);

        % p for Age^2 term
        p_val = p_age2(i,j);
        if ~isfinite(p_val) || p_val <= 0
            p_val = min_p_safe;
        end
        score = -log(p_val);
        if score>cap_score, score = cap_score; end

        bin = assign_bin(score);
        thickness_val = thickness_map(bin);

        coeff2 = beta_age2(i,j);  % quadratic coefficient
        if coeff2 >= 0
            color = quad_pos_colors(bin,:);
        else
            color = quad_neg_colors(bin,:);
        end

        rA = genesInfo(i,:); rB = genesInfo(j,:);
        fprintf(fid, '%s %d %d %s %d %d color=%d,%d,%d,thickness=%.2f\n', ...
            rA.Network{1}, rA.Start, rA.End, ...
            rB.Network{1}, rB.Start, rB.End, ...
            color(1), color(2), color(3), thickness_val);
    end
    fclose(fid);
else
    file_quad = '(no quadratic edges written)';
end

fprintf('[Circos] Link files written:\n  Linear   : %s\n  Quadratic: %s\n', file_lin, file_quad);

% (Optional) brief bin count summary
if n_lin>0
    lin_scores = -log(max(p_age(sig_lin_LT), min_p_safe));
    lin_bins = arrayfun(assign_bin, lin_scores);
    fprintf('[Circos] Linear bin counts (1..5) : %s\n', mat2str(histcounts(lin_bins,0.5:1:5.5)));
end
if n_quad>0
    quad_scores = -log(max(p_age2(sig_quad_LT), min_p_safe));
    quad_bins = arrayfun(assign_bin, quad_scores);
    fprintf('[Circos] Quadratic bin counts (1..5): %s\n', mat2str(histcounts(quad_bins,0.5:1:5.5)));
end

%% ===================== Save structured results =====================
AgeEffects = struct();
AgeEffects.Meta = struct('N_subjects',N,'alpha_fdr',alpha_fdr, ...
    'data_root',data_root,'mat_file',mat_file,'demo_csv',demo_csv, ...
    'network_index_csv',network_index_csv);
AgeEffects.Global = struct('Linear',struct('Model',mdl_lin,'AIC',AIC_lin), ...
                           'Quadratic',struct('Model',mdl_quad,'AIC',AIC_quad), ...
                           'BestModel',best_model);
AgeEffects.Edgewise = struct('beta_age',beta_age,'beta_age2',beta_age2, ...
    't_age',t_age,'t_age2',t_age2,'p_age',p_age,'p_age2',p_age2, ...
    'q_value',q_mat,'model_used',model_used, ...
    'AIC_linear',AIC_L_mat,'AIC_quadratic',AIC_Q_mat, ...
    'sig_mask_all_q',sig_mask_all,'sig_mask_linear_q',sig_mask_linear, ...
    'sig_mask_quadratic_q',sig_mask_quadratic,'peak_age_years',peak_age_yrs);

save(fullfile(out_dir,'AgeEffects_GMGM_SALD.mat'), 'AgeEffects','-v7.3');
fprintf('\nAll GM-GM results saved to: %s\n', fullfile(out_dir,'AgeEffects_GMGM_SALD.mat'));