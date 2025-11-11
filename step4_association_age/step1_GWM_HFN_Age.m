% Analysis_CorrelationWithAge_GWMHFN_latest.m
% Age effects for GWM-HFN (GM-WM-GM) in SALD with formula-free fitlm calls
% - Robust data import (Sex/Hand safe)
% - Global mean: linear vs quadratic (via explicit X matrices, no formula strings)
% - Edgewise model selection by AIC (linear vs quadratic), BH-FDR
% - Network-level proportions (within vs between) with two-proportion test + Wilson CI
% - Circos link files (linear and quadratic)
% - Structured outputs for downstream comparison with GM-GM (step2)
%
% Author: Zeqiang Linli
% Date: 2025-08-13

clear; clc; close all;

%% ===================== User configuration =====================
data_root = 'E:\Neuroimage\MyProject\GMWM_Network\Data\SALD';
mat_file  = fullfile(data_root, 'FCMData_GWMHFN.mat');          % contains fcm_GWG_cstd_prod, gmean_fcm_GWG_cstd_prod
demo_csv  = fullfile(data_root, 'Anlysis_DemographicInfo.csv'); % demographics

% Network index (AAL90 7-network assignment)
network_index_csv = 'E:\Neuroimage\MyProject\GMWM_Network\Data\AAL_7networkIndex.csv';
network_names = {'VN','SMN','AN','LB','FPN','DMN','BGN'};

% Circos genes info (AAL)
genes_info_file = 'E:\Neuroimage\Software\circos-0.69-6\AAL\genes.labelsAAL.links.txt';

% Output directory
out_dir = 'E:\Neuroimage\MyProject\GMWM_Network\Code_R1\step4_association_age';
fig_dir = fullfile('E:\Neuroimage\MyProject\GMWM_Network\Code_R1\visualization\Figure6_Age');
circos_dir = fullfile(out_dir, 'circos');
if ~exist(out_dir,'dir');   mkdir(out_dir);   end
if ~exist(fig_dir,'dir');   mkdir(fig_dir);   end
if ~exist(circos_dir,'dir');mkdir(circos_dir);end

% FDR alpha for edgewise significance and visualization
alpha_fdr = 0.05;


%% ===================== Load connectivity and demographics =====================
% Load matrices
S = load(mat_file);
if ~isfield(S, 'fcm_GWG_cstd_prod') || ~isfield(S, 'gmean_fcm_GWG_cstd_prod')
    error('Expected variables fcm_GWG_cstd_prod (90x90xN) and gmean_fcm_GWG_cstd_prod (Nx1) in %s', mat_file);
end
Z_GWMHFN        = S.fcm_GWG_cstd_prod;            % 90x90xN
Global_GWMHFN   = S.gmean_fcm_GWG_cstd_prod(:);   % Nx1
[n1, n2, N] = size(Z_GWMHFN);
if n1~=90 || n2~=90, error('Z_GWMHFN must be 90x90xN'); end
if numel(Global_GWMHFN)~=N, error('Global vector length does not match N'); end

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
AgeZ = zscore(Age);
FDZ  = zscore(FD);
AgeZ2 = AgeZ.^2; % quadratic term explicitly

%% ===================== Global mean analysis (no formula strings) =====================
% Build design matrices
Xg_lin  = [AgeZ,     SexM, HandR, FDZ];
Xg_quad = [AgeZ,AgeZ2,SexM, HandR, FDZ];
Yg = Global_GWMHFN;

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

fprintf('\n=== Global mean GWM-HFN analysis ===\n');
fprintf('AIC (Linear)   = %.3f\n', AIC_lin);
fprintf('AIC (Quadratic)= %.3f\n', AIC_quad);
fprintf('Best model     = %s\n', best_model);

%% ===================== Global mean scatter plot (Linear model preferred) =====================
% Draw scatter (Age vs Global GWM-HFN) with density, regression line, and partial correlation.
% Note:
% - Requires "scatplot" function available on MATLAB path.
% - Saves figure to fig_dir.

if strcmpi(best_model, 'Linear')
    Y = Global_GWMHFN;  % Global mean of GWM-HFN

    % Quick check for NaNs to avoid issues in plotting/regression
    valid_mask = ~isnan(Age) & ~isnan(Y);
    % Figure settings
    figure('Position', [100, 100, 800, 600]);
    set(gcf, 'Color', 'white');

    % Density scatter
    out = scatplot(Age(valid_mask), Y(valid_mask), [], [], [], [], [], 25);
    colorbar('off');
    hold on;

    % Linear fit line (same as reference: fit Age -> Y and predict on 18~80)
    mdl_sc = fitlm(Age(valid_mask), Y(valid_mask));  % simple linear model for the overlay
    x_range = linspace(18, 80, 100);
    y_fit = predict(mdl_sc, x_range');
    [y_fit, y_ci] = predict(mdl_sc, x_range');  % y_ci(:,1)=下限, y_ci(:,2)=上限
    fill([x_range fliplr(x_range)], [y_ci(:,1)' fliplr(y_ci(:,2)')], ...
        [0.8 0.2 0.2], 'FaceAlpha', 0.18, 'EdgeColor','none');
    plot(x_range, y_fit, 'LineWidth', 2.5, 'Color', [0.8 0.2 0.2]);  % red line

    % Axes limits and ticks (kept as in the reference script)
    xlim([18, 80]);
    ylim([0, 50]);
    xticks([20, 30, 40, 50, 60, 70, 80]);
    yticks(5:10:50);

    % Partial correlation controlling Sex, Hand, and FD (as in reference X(:,2:4))
    % Here we explicitly pass covariates: SexM, HandR, FD
    covars = [SexM, HandR, FD];
    % Align covariates with valid_mask
    covars = covars(valid_mask, :);
    [r, p] = partialcorr(Y(valid_mask), Age(valid_mask), covars, 'Rows', 'complete');

    % Statistical annotation
    stats_str = sprintf('r = %.3f, p = %.1e', r, p);
    text(0.05, 0.95, stats_str, 'Units', 'normalized', ...
        'FontSize', 20, 'FontWeight', 'bold', ...
        'BackgroundColor', [1 1 1 0.8], ...
        'EdgeColor', [0.5 0.5 0.5], ...
        'Margin', 3);

        % Global aesthetics and grid (as in reference)
    set(gca, 'FontSize', 16, 'FontWeight', 'bold', ...
        'LineWidth', 1.5, 'Box', 'on', ...
        'TickDir', 'out', 'TickLength', [0.01 0.01], ...
        'Color', 'white');
       % Title and labels (kept verbatim)
    title('Correlation of Mean GWM-HFN and Age', ...
        'FontSize', 25, 'FontWeight', 'bold', ...
        'Units', 'normalized', 'Position', [0.5, 1.02]);

    xlabel('Age (years)', 'FontSize', 23, 'FontWeight', 'bold');
    ylabel('Mean GWM-HFN', 'FontSize', 23, 'FontWeight', 'bold');


    grid on;
    set(gca, 'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.3, 'GridLineStyle', '-');
    set(gcf, 'Position', [100, 100, 700, 600]);

    % Save figure (optimized: use fig_dir)
    out_name = fullfile(fig_dir, 'Scatterplot_meanGWMHFN_Age');
    print([out_name '.tiff'], '-dtiff', '-r300');
    saveas(gcf, [out_name '.fig']);
    print([out_name '.svg'], '-dsvg');

    % Cleanup handle (optional)
    clear mdl_sc x_range y_fit covars stats_str out valid_mask
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

% Reuse covariates matrices for speed
Xe_lin  = [AgeZ,     SexM, HandR, FDZ];           % VarNames: {'AgeZ','SexM','HandR','FDZ','Y'}
Xe_quad = [AgeZ,AgeZ2,SexM, HandR, FDZ];          % VarNames: {'AgeZ','AgeZ2','SexM','HandR','FDZ','Y'}

fprintf('\nFitting %d edges (linear vs quadratic)...\n', n_edges);
for k = 1:n_edges
    i = idx_i(k); j = idx_j(k);
    Y = squeeze(Z_GWMHFN(i,j,:));  % Nx1

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
        % AgeZ coefficient
        c1 = mdl_Q.Coefficients.Estimate('AgeZ');
        t1 = mdl_Q.Coefficients.tStat('AgeZ');
        p1 = mdl_Q.Coefficients.pValue('AgeZ');
        % AgeZ^2 coefficient
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

        beta_age(i,j)  = c1; t_age(i,j)  = t1; p_age(i,j)  = p1;
    end
end

% Combined p for FDR
combined_p = nan(90,90);
lin_sel  = (model_used==1);
quad_sel = (model_used==2);
combined_p(lin_sel)  = p_age(lin_sel);
combined_p(quad_sel) = p_age2(quad_sel);

% Benjamini-Hochberg FDR (inline)
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

% Significant masks (q < alpha)
sig_mask_all = false(90,90);
sig_mask_all(maskLT) = q_vec < alpha_fdr;
sig_mask_all = sig_mask_all | sig_mask_all.'; % symmetric copy for convenience

sig_mask_linear    = (model_used==1) & (q_mat<alpha_fdr);
sig_mask_quadratic = (model_used==2) & (q_mat<alpha_fdr);

fprintf('Edgewise results: %d/%d significant edges at q<%.2f.\n', ...
    sum((q_vec<alpha_fdr), 'omitnan'), nnz(maskLT), alpha_fdr);
% Edgewise results: 2044/4005 significant edges at q<0.05.

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
    'VariableNames', {'Node1','Node2','Net1','Net2','ModelUsed',...
    'Beta_Age','t_Age','p_Age',...
    'Beta_Age2','t_Age2','p_Age2',...
    'q_value','PeakAge_years'});

csv_edge = fullfile(out_dir, 'Edgewise_AgeEffects_GWMHFN.csv');
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

save(fullfile(out_dir,'Matrix_Age.mat'),  'matrix_age');
save(fullfile(out_dir,'Matrix_Age2.mat'), 'matrix_age2');

% Quick visuals
f1 = figure('Color','w'); imagesc(matrix_age); axis square; colorbar;
title('Significant Age (linear) coefficients, q<0.05'); print(f1, fullfile(fig_dir,'Matrix_Age_linear.tiff'), '-dtiff','-r300');
f2 = figure('Color','w'); imagesc(matrix_age2); axis square; colorbar;
title('Significant Age^2 (quadratic) coefficients, q<0.05'); print(f2, fullfile(fig_dir,'Matrix_Age_quadratic.tiff'), '-dtiff','-r300');

%% ===================== Network-level proportions (within vs between) =====================
% -------- 1. Basic checks and preparation --------
if ~exist('net_index','var')
    error('net_index not found in workspace.');
end
net_index = net_index(:)';          % ensure row
num_nodes = numel(net_index);
if num_nodes ~= size(model_used,1)
    error('Size mismatch: net_index length (%d) vs model_used size (%d).', ...
        num_nodes, size(model_used,1));
end

K = numel(unique(net_index));
% Remap network labels to 1..K if they are not already
[~,~,net_index] = unique(net_index);  % now exactly 1..K

% Build lower-triangle mask
maskLT = tril(true(num_nodes), -1);

% -------- 2. Build combined significant edge mask (lower triangle) --------
sig_lin_LT  = (model_used==1) & (q_mat < alpha_fdr) & maskLT;
sig_quad_LT = (model_used==2) & (q_mat < alpha_fdr) & maskLT;
sig_combined_LT = (sig_lin_LT | sig_quad_LT);  % merged significant edges

% -------- 3. Global within vs between (for reference / printout) --------
% Node counts per network
n_k = zeros(K,1);
for k = 1:K
    n_k(k) = sum(net_index==k);
end

% Build a within mask (LT) across all networks
within_mask_LT = false(num_nodes);
for k = 1:K
    nodes_k = find(net_index==k);
    tmp = false(num_nodes); tmp(nodes_k, nodes_k) = true;
    within_mask_LT = within_mask_LT | (tmp & maskLT);
end
between_mask_LT = maskLT & ~within_mask_LT;

% Global denominators
Total_edges          = num_nodes*(num_nodes-1)/2;         % C(N,2)
Den_within_global    = sum(n_k.*(n_k-1)/2);               % Σ C(n_k,2)
Den_between_global   = Total_edges - Den_within_global;

% Global numerators
Num_within_global    = nnz(sig_combined_LT & within_mask_LT);
Num_between_global   = nnz(sig_combined_LT & between_mask_LT);

Prop_within_global   = Num_within_global  / max(Den_within_global, 1);
Prop_between_global  = Num_between_global / max(Den_between_global,1);

fprintf('\n=== Combined (Linear+Quadratic) global proportions (q < %.3f) ===\n', alpha_fdr);
fprintf('Within-network:  %d / %d = %.4f\n', Num_within_global, Den_within_global, Prop_within_global);
fprintf('Between-network: %d / %d = %.4f\n', Num_between_global, Den_between_global, Prop_between_global);

% -------- 4. Per-network numerators (within & between) --------
% Extract edge list (lower triangle)
[i_edges, j_edges] = find(sig_combined_LT);     % i_edges > j_edges (by construction of maskLT)
net_i = net_index(i_edges);
net_j = net_index(j_edges);

within_count  = zeros(K,1);
between_count = zeros(K,1);

% For each significant edge, update counts
for e = 1:numel(i_edges)
    if net_i(e) == net_j(e)
        % Within-network edge contributes to exactly one network
        within_count(net_i(e)) = within_count(net_i(e)) + 1;
    else
        % Between-network edge contributes to BOTH endpoint networks
        between_count(net_i(e)) = between_count(net_i(e)) + 1;
        between_count(net_j(e)) = between_count(net_j(e)) + 1;
    end
end

% -------- 5. Per-network denominators --------
den_within_k  = n_k.*(n_k-1)/2;           % C(n_k,2)
den_between_k = n_k .* (num_nodes - n_k); % n_k*(N - n_k)

% Proportions
prop_within_k  = within_count  ./ max(den_within_k, 1);
prop_between_k = between_count ./ max(den_between_k,1);

% Replace NaNs if any (e.g., n_k=1 -> den_within=0)
prop_within_k(~isfinite(prop_within_k))   = 0;
prop_between_k(~isfinite(prop_between_k)) = 0;

% -------- 6. Figure: per-network grouped bar (Within vs Between) --------
if ~exist('network_names','var') || numel(network_names)~=K
    network_names = arrayfun(@(x) sprintf('Net%d', x), 1:K, 'UniformOutput', false);
end

figure('Position',[120 120 750 650]); set(gcf,'Color','white');

bar_data = [prop_within_k, prop_between_k];  % K x 2
hB = bar(bar_data, 'grouped');
hold on;

% Optional color customization (matching earlier palette, e.g., blue vs orange)
if size(hB,2) == 2
    set(hB(1), 'FaceColor', [76 114 176]/255);   % Within
    set(hB(2), 'FaceColor', [221 132 82]/255);   % Between
end

set(gca, 'XTickLabel', network_names, ...
    'FontName','Arial','FontSize',16,'FontWeight','bold',...
    'LineWidth',1.5,'TickDir','out','Box','on');
ylabel('Proportion of Significant Edges', 'FontSize',22,'FontWeight','bold');
title('Age-Related Edges by Networks', 'FontSize',24,'FontWeight','bold');

legend({'Within-Network','Between-Network'},'Location','northwest','FontSize',13,'Box','off');

% Annotate percentages above bars
numGroups = K;
numBars   = 2;
groupWidth = min(0.8, numBars/(numBars+1.5));

for b = 1:numBars
    % X positions of the bars in group b
    xCenters = (1:numGroups) - groupWidth/2 + (2*b-1)*groupWidth/(2*numBars);
    for k = 1:numGroups
        val = bar_data(k,b);
        text(xCenters(k), val, sprintf('%.1f%%', val*100), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',14,'FontWeight','bold');
    end
end

grid on;
set(gca,'YGrid','on','GridAlpha',0.15,'GridColor',[0.7 0.7 0.7]);

% (Optional) y-limit for nicer spacing
ylim([0, max(bar_data(:))*1.15]);

% Save
out_file = fullfile(fig_dir,'PerNetwork_WithinBetween_Combined.tiff');
print(gcf, out_file, '-dtiff','-r300');
out_file = fullfile(fig_dir,'PerNetwork_WithinBetween_Combined.svg');
print(gcf, out_file, '-dsvg');

fprintf('\nPer-network summary (first 5 rows shown):\n');
perNet_table = table(network_names(:), n_k, within_count, den_within_k, prop_within_k, ...
                                  between_count, den_between_k, prop_between_k, ...
    'VariableNames', {'Network','Nnodes','WithinNum','WithinDen','WithinProp','BetweenNum','BetweenDen','BetweenProp'});
writetable(perNet_table, fullfile(fig_dir,'PerNetwork_WithinBetween_Combined.xlsx'));

disp(table(network_names(:), n_k, within_count, between_count, ...
    prop_within_k, prop_between_k, 'VariableNames', ...
    {'Network','N_nodes','Within_Num','Between_Num','Within_Prop','Between_Prop'}));
%% ===================== Circos link files =====================
% Read genes info (tab-delimited: Network, Start, End, Region)
genesInfo = readtable(genes_info_file, 'FileType','text', 'Delimiter','\t', 'ReadVariableNames', false);
genesInfo.Properties.VariableNames = {'Network','Start','End','Region'};
if height(genesInfo)<90
    error('genesInfo has less than 90 entries');
end

% Prepare vectors
LT_idx   = find(maskLT);
q_valsLT = q_mat(maskLT);
neglog10q = -log10(q_valsLT);

% Linear links
sel_lin = sig_mask_linear(maskLT);
rows_lin = idx_i(sel_lin);
cols_lin = idx_j(sel_lin);
coeff_lin = beta_age(maskLT); coeff_lin = coeff_lin(sel_lin);
nlq_lin = neglog10q(sel_lin);

link_file_lin = fullfile(circos_dir, 'segdup_GWMHFN_age_linear.link.txt');
fid = fopen(link_file_lin, 'w');
for k = 1:numel(rows_lin)
    i = rows_lin(k); j = cols_lin(k);
    rA = genesInfo(i,:); rB = genesInfo(j,:);
    c  = coeff_lin(k);
    if c>0, color = [255, 80, 80]; else, color = [55, 140, 231]; end
    nlq = nlq_lin(k);
    if isnan(nlq) || nlq<=1, th=1;
    elseif nlq<=2, th=2;
    elseif nlq<=3, th=3;
    elseif nlq<=4, th=5;
    else, th=8; end
    fprintf(fid, '%s %d %d %s %d %d color=%d,%d,%d,thickness=%.2f\n', ...
        rA.Network{1}, rA.Start, rA.End, ...
        rB.Network{1}, rB.Start, rB.End, ...
        color(1), color(2), color(3), th);
end
fclose(fid);

% Quadratic links
sel_quad = sig_mask_quadratic(maskLT);
rows_quad = idx_i(sel_quad);
cols_quad = idx_j(sel_quad);
coeff_quad = beta_age2(maskLT); coeff_quad = coeff_quad(sel_quad);
nlq_quad = neglog10q(sel_quad);

link_file_quad = fullfile(circos_dir, 'segdup_GWMHFN_age_quadratic.link.txt');
fid = fopen(link_file_quad, 'w');
for k = 1:numel(rows_quad)
    i = rows_quad(k); j = cols_quad(k);
    rA = genesInfo(i,:); rB = genesInfo(j,:);
    c  = coeff_quad(k);
    if c>0, color = [162, 148, 249]; else, color = [128, 157, 60]; end
    nlq = nlq_quad(k);
    if isnan(nlq) || nlq<=1, th=1;
    elseif nlq<=2, th=2;
    elseif nlq<=3, th=3;
    elseif nlq<=4, th=5;
    else, th=8; end
    fprintf(fid, '%s %d %d %s %d %d color=%d,%d,%d,thickness=%.2f\n', ...
        rA.Network{1}, rA.Start, rA.End, ...
        rB.Network{1}, rB.Start, rB.End, ...
        color(1), color(2), color(3), th);
end
fclose(fid);

fprintf('\nCircos link files saved:\n  %s\n  %s\n', link_file_lin, link_file_quad);

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


save(fullfile(out_dir,'AgeEffects_GWMHFN_SALD.mat'), 'AgeEffects','-v7.3');
fprintf('\nAll results saved to: %s\n', fullfile(out_dir,'AgeEffects_GWMHFN_SALD.mat'));