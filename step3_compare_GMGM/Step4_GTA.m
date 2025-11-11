% MATLAB Code for Comparing Network Properties (GWM-HFN vs GM-GM)
% --- VERSION 3: Robust SE, Clean Visuals, FDR on Global AUC, and Full SourceData Exports ---

clear; clc; close all;


%% --- Configuration ---
% 1. File Paths and Naming
path_gwmhfn = 'E:\Neuroimage\MyProject\GMWM_Network\Data\SLIM\NetworkAnalysis\NetworkResults\time1';
path_gmgm   = 'E:\Neuroimage\MyProject\GMWM_Network\Data\SLIM\\NetworkAnalysis_GMGM\NetworkResults\time1';
file_prefix_gwmhfn = 'GTA_GMWMGM_';
file_prefix_gmgm   = 'GTA_GMGM_';
file_suffix = '.mat';

% 2. Network Property Names
name_global = {'CpSet','LpSet','GammaSet','LambdaSet','SigmaSet','ElocSet','EglobSet','AssortativitySet','ModularitySet'};
name_local  = {'DegreeSet','NodalEfficiencySet','BetweennessSet','ClusteringCoefficientSet','ParticipantCoefficientSet','SubgraphCentralitySet','EigenvectorCentralitySet','PageRankCentralitySet'};

num_global_metrics = length(name_global);
num_local_metrics  = length(name_local);
num_nodes = 90;

% 3. Sparsity Settings
sparsity_vector = 0.01:0.01:0.50;
num_sparsities  = length(sparsity_vector);
auc_sparsities_range  = [0.10, 0.34];
auc_indices           = find(sparsity_vector >= auc_sparsities_range(1) & sparsity_vector <= auc_sparsities_range(2));
auc_sparsities_values = sparsity_vector(auc_indices);

% 4. Statistical Settings
alpha_level = 0.05;             % global default
alpha_level_nodal_viz = 0.001;  % for Figure 4 visualization threshold (kept as in your code)

% 5. Figure/Color defaults
color_gwmhfn = [0, 0.4470, 0.7410]; % Blue
color_gmgm   = [0.8500, 0.3250, 0.0980]; % Red
color_sig    = [0.2, 0.2, 0.2]; % Dark gray
color_text   = [0.2, 0.2, 0.2]; % Dark gray

% 6. Output directories
output_dir = fullfile(pwd, 'Results_GTA');
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

% publication figure dir 
fig_dir = 'E:\Neuroimage\MyProject\GMWM_Network\Code_R1\visualization\Figure5_CompareGMGM';
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

% 7. SourceData Excel (保存所有绘图的原始数据)
excel_file = fullfile(fig_dir, 'Figure5_vsGMGM_SourceData.xlsx'); % 您可改到固定目录

% Helper to create sheet names <= 31 chars
mkSheet = @(s) s(1:min(31, numel(s)));

%% --- Get Subject List ---
gwmhfn_files = dir(fullfile(path_gwmhfn, [file_prefix_gwmhfn, '*', file_suffix]));
subject_list = cell(length(gwmhfn_files), 1);
for i = 1:length(gwmhfn_files)
    fname = gwmhfn_files(i).name;
    subname_part = fname(length(file_prefix_gwmhfn)+1 : end-length(file_suffix));
    subject_list{i} = subname_part;
end
num_subjects = length(subject_list);
fprintf('Found %d subjects.\n', num_subjects);

%% --- Data Loading ---
fprintf('Loading data...\n');
all_data_gwmhfn.global = struct(); all_data_gwmhfn.local = struct();
all_data_gmgm.global   = struct(); all_data_gmgm.local   = struct();

for i = 1:num_global_metrics
    all_data_gwmhfn.global.(name_global{i}) = NaN(num_subjects, num_sparsities);
    all_data_gmgm.global.(name_global{i})   = NaN(num_subjects, num_sparsities);
end
for i = 1:num_local_metrics
    all_data_gwmhfn.local.(name_local{i}) = NaN(num_subjects, num_sparsities, num_nodes);
    all_data_gmgm.local.(name_local{i})   = NaN(num_subjects, num_sparsities, num_nodes);
end

for s = 1:num_subjects
    sub_id = subject_list{s};
    fname_gwmhfn = fullfile(path_gwmhfn, [file_prefix_gwmhfn, sub_id, file_suffix]);
    fname_gmgm   = fullfile(path_gmgm,   [file_prefix_gmgm,   sub_id, file_suffix]);
    if exist(fname_gwmhfn, 'file') && exist(fname_gmgm, 'file')
        data_gwmhfn = load(fname_gwmhfn);
        data_gmgm   = load(fname_gmgm);
        % Global metrics
        for i = 1:num_global_metrics
            metric_name = name_global{i};
            if isfield(data_gwmhfn, metric_name) && isfield(data_gmgm, metric_name)
                all_data_gwmhfn.global.(metric_name)(s, :) = data_gwmhfn.(metric_name);
                all_data_gmgm.global.(metric_name)(s, :)   = data_gmgm.(metric_name);
            else
                warning('Metric %s not found for subject %s', metric_name, sub_id);
            end
        end
        % Local metrics
        for i = 1:num_local_metrics
            metric_name = name_local{i};
            if isfield(data_gwmhfn, metric_name) && isfield(data_gmgm, metric_name)
                gwmhfn_local_data = squeeze(data_gwmhfn.(metric_name));
                gmgm_local_data   = squeeze(data_gmgm.(metric_name));
                % GWM-HFN orientation
                if size(gwmhfn_local_data, 2) == num_nodes && size(gwmhfn_local_data, 1) == num_sparsities
                    gwmhfn_local_data = gwmhfn_local_data';
                elseif ~(size(gwmhfn_local_data, 1) == num_nodes && size(gwmhfn_local_data, 2) == num_sparsities)
                    warning('Unexpected local data dimensions for %s subject %s (GWM-HFN)', metric_name, sub_id);
                    continue;
                end
                all_data_gwmhfn.local.(metric_name)(s, :, :) = gwmhfn_local_data';
                % GM-GM orientation
                if size(gmgm_local_data, 2) == num_nodes && size(gmgm_local_data, 1) == num_sparsities
                    gmgm_local_data = gmgm_local_data';
                elseif ~(size(gmgm_local_data, 1) == num_nodes && size(gmgm_local_data, 2) == num_sparsities)
                    warning('Unexpected local data dimensions for %s subject %s (GM-GM)', metric_name, sub_id);
                    continue;
                end
                all_data_gmgm.local.(metric_name)(s, :, :) = gmgm_local_data';
            else
                warning('Metric %s not found for subject %s', metric_name, sub_id);
            end
        end
    else
        warning('Data files not found for subject %s', sub_id);
    end
end
fprintf('Data loading complete.\n');

%% --- Global Property Analysis ---
fprintf('\n--- Analyzing Global Properties ---\n');
results_global = struct();
all_global_auc_pvals = NaN(num_global_metrics, 1); % Collect p-values for FDR across metrics

for i = 1:num_global_metrics
    metric_name = name_global{i};
    fprintf('Analyzing %s...\n', metric_name);

    data1 = all_data_gwmhfn.global.(metric_name); % [subjects x sparsities]
    data2 = all_data_gmgm.global.(metric_name);   % [subjects x sparsities]

    % 1) Point-wise paired t-tests (uncorrected) for 0.10–0.34 (for curve asterisks)
    p_vals_pointwise = NaN(length(auc_indices), 1);
    t_stats_pointwise = NaN(length(auc_indices), 1);
    for k = 1:length(auc_indices)
        idx = auc_indices(k);
        diff_data = data1(:, idx) - data2(:, idx);
        valid = ~isnan(diff_data);
        if sum(valid) > 1
            [~, p, ~, stats] = ttest(diff_data(valid));
            p_vals_pointwise(k) = p;
            t_stats_pointwise(k) = stats.tstat;
        end
    end
    results_global.(metric_name).pointwise_p_uncorrected = p_vals_pointwise;
    results_global.(metric_name).pointwise_t = t_stats_pointwise;
    results_global.(metric_name).pointwise_sparsities = auc_sparsities_values;

    % 2) AUC per subject (0.10–0.34)
    auc_data1 = NaN(num_subjects, 1);
    auc_data2 = NaN(num_subjects, 1);
    for s = 1:num_subjects
        v1 = data1(s, auc_indices); v2 = data2(s, auc_indices);
        valid = ~isnan(v1) & ~isnan(v2);
        if sum(valid) > 1
            auc_data1(s) = trapz(auc_sparsities_values(valid), v1(valid));
            auc_data2(s) = trapz(auc_sparsities_values(valid), v2(valid));
        end
    end
    results_global.(metric_name).auc_gwmhfn_all_subjects = auc_data1;
    results_global.(metric_name).auc_gmgm_all_subjects   = auc_data2;

    % 3) AUC paired t-test
    diff_auc = auc_data1 - auc_data2;
    valid_auc = ~isnan(diff_auc);
    if sum(valid_auc) > 1
        [h_auc, p_auc, ci_auc, stats_auc] = ttest(diff_auc(valid_auc));
        cohen_d_auc = cohensD_paired(auc_data1(valid_auc), auc_data2(valid_auc));

        all_global_auc_pvals(i) = p_auc;
        results_global.(metric_name).auc_p_uncorrected = p_auc;
        results_global.(metric_name).auc_t = stats_auc.tstat;
        results_global.(metric_name).auc_df = stats_auc.df;
        results_global.(metric_name).auc_cohen_d = cohen_d_auc;
        results_global.(metric_name).auc_significant_uncorrected = h_auc;
        results_global.(metric_name).auc_mean_gwmhfn = mean(auc_data1, 'omitnan');
        results_global.(metric_name).auc_se_gwmhfn   = std(auc_data1, 'omitnan')/sqrt(sum(~isnan(auc_data1)));
        results_global.(metric_name).auc_mean_gmgm   = mean(auc_data2, 'omitnan');
        results_global.(metric_name).auc_se_gmgm     = std(auc_data2, 'omitnan')/sqrt(sum(~isnan(auc_data2)));

        fprintf('  AUC comparison (uncorrected): p = %.4g, t(%d) = %.3f, d = %.3f\n', p_auc, stats_auc.df, stats_auc.tstat, cohen_d_auc);
    else
        fprintf('  AUC comparison: Not enough valid pairs.\n');
        results_global.(metric_name).auc_p_uncorrected = NaN;
        results_global.(metric_name).auc_significant_uncorrected = NaN;
    end
end

% 4) FDR correction across global metrics (AUC p-values)
valid_p_idx = ~isnan(all_global_auc_pvals);
if any(valid_p_idx)
    q_vals_global_auc_fdr = mafdr(all_global_auc_pvals(valid_p_idx), 'BHFDR', true);
    idxs = find(valid_p_idx);
    for k = 1:length(idxs)
        metric_name = name_global{idxs(k)};
        results_global.(metric_name).auc_q_fdr_across_metrics = q_vals_global_auc_fdr(k);
        results_global.(metric_name).auc_significant_fdr = q_vals_global_auc_fdr(k) < alpha_level;
        fprintf('  FDR(AUC) %s: q = %.4g (sig: %d)\n', metric_name, q_vals_global_auc_fdr(k), q_vals_global_auc_fdr(k) < alpha_level);
    end
else
    fprintf('No valid AUC p-values for FDR.\n');
end
  % FDR(AUC) CpSet: q = 3.251e-107 (sig: 1)
  % FDR(AUC) LpSet: q = 2.112e-133 (sig: 1)
  % FDR(AUC) GammaSet: q = 0.2114 (sig: 0)
  % FDR(AUC) LambdaSet: q = 9.632e-136 (sig: 1)
  % FDR(AUC) SigmaSet: q = 1.813e-37 (sig: 1)
  % FDR(AUC) ElocSet: q = 1.743e-07 (sig: 1)
  % FDR(AUC) EglobSet: q = 2.922e-159 (sig: 1)
  % FDR(AUC) AssortativitySet: q = 2.682e-100 (sig: 1)
  % FDR(AUC) ModularitySet: q = 0.3406 (sig: 0)
%% --- Local Property Analysis (AUC per node with FDR across nodes) ---
fprintf('\n--- Analyzing Local Properties ---\n');
results_local = struct();
for i = 1:num_local_metrics
    metric_name = name_local{i};
    fprintf('Analyzing %s...\n', metric_name);

    data1 = all_data_gwmhfn.local.(metric_name); % [subs x sparsity x nodes]
    data2 = all_data_gmgm.local.(metric_name);

    p_vals_nodal_auc      = NaN(num_nodes, 1);
    t_stats_nodal_auc     = NaN(num_nodes, 1);
    cohen_d_nodal_auc     = NaN(num_nodes, 1);
    mean_auc_gwmhfn_nodal = NaN(num_nodes, 1);
    mean_auc_gmgm_nodal   = NaN(num_nodes, 1);

    for n = 1:num_nodes
        node1 = squeeze(data1(:, auc_indices, n));
        node2 = squeeze(data2(:, auc_indices, n));
        auc1 = NaN(num_subjects, 1);
        auc2 = NaN(num_subjects, 1);
        for s = 1:num_subjects
            v1 = node1(s, :); v2 = node2(s, :);
            valid = ~isnan(v1) & ~isnan(v2);
            if sum(valid) > 1
                auc1(s) = trapz(auc_sparsities_values(valid), v1(valid));
                auc2(s) = trapz(auc_sparsities_values(valid), v2(valid));
            end
        end
        diff_auc = auc1 - auc2;
        valid = ~isnan(diff_auc);
        if sum(valid) > 1
            [~, p, ~, stats] = ttest(diff_auc(valid));
            p_vals_nodal_auc(n)  = p;
            t_stats_nodal_auc(n) = stats.tstat;
            cohen_d_nodal_auc(n) = cohensD_paired(auc1(valid), auc2(valid));
            mean_auc_gwmhfn_nodal(n) = mean(auc1, 'omitnan');
            mean_auc_gmgm_nodal(n)   = mean(auc2, 'omitnan');
        end
    end
    q_vals_nodal_auc_fdr = NaN(size(p_vals_nodal_auc));
    valid = ~isnan(p_vals_nodal_auc);
    if any(valid)
        q_vals_nodal_auc_fdr(valid) = mafdr(p_vals_nodal_auc(valid), 'BHFDR', true);
    end

    results_local.(metric_name).nodal_auc_p   = p_vals_nodal_auc;
    results_local.(metric_name).nodal_auc_t   = t_stats_nodal_auc;
    results_local.(metric_name).nodal_auc_q_fdr = q_vals_nodal_auc_fdr;
    results_local.(metric_name).nodal_auc_cohen_d = cohen_d_nodal_auc;
    results_local.(metric_name).nodal_mean_auc_gwmhfn = mean_auc_gwmhfn_nodal;
    results_local.(metric_name).nodal_mean_auc_gmgm   = mean_auc_gmgm_nodal;

    fprintf('Done (FDR across nodes).\n');
end

%% --- Figure: Selected global metrics boxplot  ---
fprintf('\n--- Visualizing Selected Global Metrics (Boxplot) ---\n');
figure('Name', 'Selected Global Metrics Boxplot', 'Position', [50, 50, 900, 300], 'Color', 'white');

% Define the metrics to plot in specified order
selected_metrics = {'CpSet', 'AssortativitySet', 'EglobSet', 'LpSet', 'LambdaSet', 'SigmaSet'};
num_selected = length(selected_metrics);

% Define colors
color_gwmhfn = [0.8500, 0.3250, 0.0980]; % Red
color_gmgm = [0, 0.4470, 0.7410]; % Blue
color_text = [0.2, 0.2, 0.2]; % Dark gray for text

ax_handles = gobjects(1, num_selected); % Preallocate array for axes handles

% Create 1x6 subplot layout
for i = 1:num_selected
    metric_name = selected_metrics{i};
    ax_handles(i) = subplot(1, num_selected, i);

    subplot(1, num_selected, i);

    % Prepare data for boxplot
    data1 = all_data_gwmhfn.global.(metric_name)(:, auc_indices);
    data2 = all_data_gmgm.global.(metric_name)(:, auc_indices);

    % Flatten AUC data
    plot_data = [data1(:); data2(:)];
    plot_group = [ones(numel(data1), 1); 2*ones(numel(data2), 1)];

    % Create boxplot
    boxplot(plot_data, plot_group, 'Colors', [color_gwmhfn; color_gmgm], 'Symbol', 'k.');

    % Add significance information if available
    if isfield(results_global.(metric_name), 'auc_p_uncorrected') && ...
            ~isnan(results_global.(metric_name).auc_p_uncorrected)
        % p_val_text = sprintf('p = %.3e', results_global.(metric_name).auc_p_uncorrected);
        % d_val_text = sprintf('d = %.2f', results_global.(metric_name).auc_cohen_d);
        q_val_text = '';
        if isfield(results_global.(metric_name), 'auc_q_fdr_across_metrics') && ~isnan(results_global.(metric_name).auc_q_fdr_across_metrics)
            q_val_text = sprintf('p_{FDR} = %.3e', results_global.(metric_name).auc_q_fdr_across_metrics);
        end
        ylim_curr = ylim;
        text(1.5, ylim_curr(2)*0.98, sprintf('%s', q_val_text), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 8,'FontWeight', 'bold', ...
            'Color', color_text, 'BackgroundColor', 'white', 'Margin', 1);
    end

    % Customize plot
    set(gca, 'XTickLabel', {'GWM-HFN', 'GM-GM'},'FontSize', 7,'FontWeight', 'bold');
    ylabel('AUC Value', 'FontSize', 10, 'FontWeight', 'bold');
    title(strrep(metric_name, 'Set', ''), 'Interpreter', 'none', 'FontSize', 12, 'FontWeight', 'bold');
    box off;
    grid on;
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);
end

sgtitle('Comparison of Global Properties (Sparsity 0.10-0.34 Integrated)', ...
    'FontSize', 14, 'FontWeight', 'bold', 'Color', color_text, 'Interpreter', 'none');

% --- Parameters to Adjust ---
% Vertical adjustments (for top margin/title)
vertical_shift = 0.05; % Shift subplots down (increase for more top space)
height_scale = 0.95;  % Scale subplot height (decrease for more top space, makes plots shorter)

% Horizontal adjustments (for side margins and inter-subplot gaps)
left_margin = 0.05;   % Space on the left of the first subplot
right_margin = 0.03;  % Space on the right of the last subplot
horizontal_gap = 0.04; % Space between adjacent subplots (increase for wider gaps)
% --- End Parameters ---

% Calculate available width and width per subplot
total_horizontal_space = 1.0 - left_margin - right_margin;
total_gap_width = (num_selected - 1) * horizontal_gap;
subplot_width = (total_horizontal_space - total_gap_width) / num_selected;

% Check if calculated width is positive
if subplot_width <= 0
    error('Calculated subplot width is negative or zero. Reduce margins or gaps, or increase figure width.');
end

% Apply adjustments to each subplot
original_bottoms = zeros(1, num_selected); % Store original bottom positions if needed later
original_heights = zeros(1, num_selected); % Store original heights

for i = 1:num_selected
    ax = ax_handles(i);
    pos = get(ax, 'Position'); % Get original position
    
    original_bottoms(i) = pos(2);
    original_heights(i) = pos(4);
    
    % Calculate new vertical position (bottom and height)
    new_height = original_heights(i) * height_scale; % Use original height for scaling
    new_bottom = original_bottoms(i) - vertical_shift; % Shift down from original position
    
    % Calculate new horizontal position (left)
    new_left = left_margin + (i - 1) * (subplot_width + horizontal_gap);
    
    % Set the new position [left, bottom, width, height]
    set(ax, 'Position', [new_left, new_bottom, subplot_width, new_height]);
end
% --- MODIFICATION END ---

fprintf('Selected global metrics boxplots saved in multiple formats to: %s\n', output_dir);

% print(fullfile(output_dir, 'Selected_Global_AUC_Boxplots.tiff'), '-dtiff', '-r600');
% print(fullfile(output_dir, 'Selected_Global_AUC_Boxplots.pdf'),  '-dpdf',  '-r600');
print(fullfile(fig_dir, 'Selected_Global_AUC_Boxplots.png'),  '-dpng',  '-r600');

fprintf('Selected global metrics boxplots rendered.\n');


%% --- Nodal Results: Save tables and visualize (Figure S2) ---
fprintf('\n--- Saving Nodal AUC Results and Visualizing ---\n');

% Load node information
node_info_file = 'E:\Neuroimage\MyProject\GMWM_Network\Data\aal_node_info.csv';
opts = detectImportOptions(node_info_file);
try
    node_info = readtable(node_info_file, opts);
    if ~isequal(sort(node_info.Node_Index), (1:90)')
        error('Node_Index in aal_node_info.csv is not 1..90.');
    end
    node_info = sortrows(node_info, 'Node_Index'); % ensure order 1..90
catch ME
    error('Failed to load node_info: %s', ME.message);
end

% Save data
for i = 1:num_local_metrics
    metric_name = name_local{i};
    q_vals = results_local.(metric_name).nodal_auc_q_fdr;
    t_vals = results_local.(metric_name).nodal_auc_t;

    Direction = strings(num_nodes,1);
    Direction(q_vals < alpha_level_nodal_viz & t_vals > 0)  = "GWM-HFN > GM-GM";
    Direction(q_vals < alpha_level_nodal_viz & t_vals < 0)  = "GM-GM > GWM-HFN";
    Direction(q_vals >= alpha_level_nodal_viz | isnan(q_vals)) = "Non-significant";

    Tnod = table((1:num_nodes)', node_info.RegionName, node_info.Network, node_info.Lobe, t_vals, q_vals, Direction, ...
        'VariableNames', {'NodeIndex','RegionName','Network','Lobe','t_value','q_value_FDR','Direction'});
    try
        writetable(Tnod, excel_file, 'Sheet', mkSheet(['FigTopo_Nodal_' metric_name]));
    catch
        writetable(Tnod, excel_file, 'Sheet', mkSheet(['Nodal_' num2str(i)]));
    end
end

% FigureS2 : Heatmap of sign-coded significance across nodes and metrics
figure('Name', 'Nodal Results Comparison', 'Position', [50, 50, 1900, 800], 'Color', 'white');

% Define visual colors
color_stronger = [0.8500, 0.3250, 0.0980]; % Red (GWM-HFN > GM-GM)
color_weaker   = [0, 0.4470, 0.7410];      % Blue (GWM-HFN < GM-GM)
color_white    = [1, 1, 1];
color_grid     = [0.8, 0.8, 0.8];

% Network labels and colors (consistent with your preference)
Networks = {'VN', 'SMN', 'AN', 'LN', 'FPN', 'DMN', 'BGN'};
network_colors = [
    186/255,  85/255, 211/255;  % VN
     70/255, 130/255, 180/255;  % SMN
    144/255, 238/255, 144/255;  % AN
    200/255,   0/255,   0/255;  % LN
    255/255, 140/255,   0/255;  % FPN
    255/255, 182/255, 193/255;  % DMN
    139/255,  69/255,  19/255   % BGN
];

% Sort nodes by network order for axis label coloring
[~, network_order] = ismember(node_info.Network, Networks);
[~, sort_idx] = sort(network_order);
node_info_sorted = node_info(sort_idx, :);

% Build results_matrix [metrics x nodes_sorted]
num_metrics = length(name_local);
results_matrix = zeros(num_metrics, num_nodes); % -1, 0, +1
for i = 1:num_metrics
    metric_name = name_local{i};
    q_vals = results_local.(metric_name).nodal_auc_q_fdr;
    t_vals = results_local.(metric_name).nodal_auc_t;
    sign_vec = zeros(1, num_nodes);
    sig_mask = (q_vals < alpha_level_nodal_viz);
    sign_vec(sig_mask & (t_vals > 0)) =  1;
    sign_vec(sig_mask & (t_vals < 0)) = -1;
    % reorder by sort_idx
    results_matrix(i, :) = sign_vec(sort_idx);
end

% Save Data
% Tm = array2table(results_matrix, 'VariableNames', cellstr(string(node_info_sorted.RegionName)), 'RowNames', strrep(name_local,'Set',''));
% try
%     writetable(Tm, excel_file, 'Sheet', mkSheet('FigTopo_NodalHeat_Sign'), 'WriteRowNames', true);
% catch
%     writetable(Tm, excel_file, 'Sheet', mkSheet('NodalHeat_Sign'), 'WriteRowNames', true);
% end

imagesc(results_matrix);
colormap([color_weaker; color_white; color_stronger]); caxis([-1 1]); hold on;
for i = 1:num_nodes+1
    plot([i-0.5, i-0.5], [0.5, num_metrics+0.5], 'Color', color_grid, 'LineWidth', 0.5);
end
for i = 1:num_metrics+1
    plot([0.5, num_nodes+0.5], [i-0.5, i-0.5], 'Color', color_grid, 'LineWidth', 0.5);
end
set(gca, 'XTick', 1:num_nodes, 'YTick', 1:num_metrics, ...
         'YTickLabel', strrep(name_local, 'Set', ''), 'FontSize', 12, ...
         'TickLength', [0 0], 'Box', 'off', 'XTickLabel', []);

% Move axes up for room for labels
set(gca, 'Position', [0.1, 0.3, 0.85, 0.6]);

% Add network separators (thicker lines)
network_boundaries = find(diff([0; findgroups(node_info.Network)]) ~= 0);
for i = 1:length(network_boundaries)
    plot([network_boundaries(i)-0.5, network_boundaries(i)-0.5], [0.5, num_metrics+0.5], ...
        'k-', 'LineWidth', 2);
end

% Bottom region labels colored by network
for i = 1:num_nodes
    net_i = node_info_sorted.Network{i};
    net_idx = find(strcmp(Networks, net_i), 1); 
    if isempty(net_idx); node_col = [0 0 0]; else; node_col = network_colors(net_idx,:); end
    text(i, num_metrics+0.55, node_info_sorted.RegionName{i}, 'Rotation', 90, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'Color', node_col, 'FontWeight', 'bold', 'FontSize', 12);
end

% Bottom network labels
for i = 1:length(Networks)
    idxs = find(strcmp(node_info_sorted.Network, Networks{i}));
    if isempty(idxs), continue; end
    midp = (min(idxs)+max(idxs))/2;
    text(midp, num_metrics+2.3, Networks{i}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'Color', network_colors(i,:), ...
        'FontWeight', 'bold', 'FontSize', 20);
end
plot([0.5, num_nodes+0.5], [num_metrics+0.5, num_metrics+0.5], 'k-', 'LineWidth', 1.5); % separator

c = colorbar('Ticks', [-1, 0, 1], 'TickLabels', {'GWM-HFN < GM-GM', 'Non-significant', 'GWM-HFN > GM-GM'});
c.Label.String = 'Direction of Difference'; c.Label.FontSize = 14; c.Label.FontWeight = 'bold';
title('Nodal Metric Comparisons', 'FontSize', 20, 'FontWeight', 'bold', 'Color', color_text);
ylabel('Local Network Metrics', 'FontSize', 18, 'FontWeight', 'bold');
plot([0.5, num_nodes+0.5], [0.5, 0.5], 'k-', 'LineWidth', 2);  % top line
plot([num_nodes+0.5, num_nodes+0.5], [0.5, num_metrics+0.5], 'k-', 'LineWidth', 2); % right line

print(fullfile(fig_dir, 'Nodal_Results_Comparison.tiff'), '-dtiff', '-r600');
print(fullfile(fig_dir, 'Nodal_Results_Comparison.pdf'),  '-dpdf',  '-r600');
print(fullfile(fig_dir, 'Nodal_Results_Comparison.png'),  '-dpng',  '-r600');
fprintf('Figure S2 saved to: %s\n', fig_dir);

fprintf('\n*** Analysis Complete. SourceData saved to: %s ***\n', excel_file);


