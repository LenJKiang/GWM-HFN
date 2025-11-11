% This script identifies hub nodes in GWM-HFN based on degree and diversity,
% and visualizes the degree distribution with fits to power law, exponential,
% and truncated power law models. It also saves the underlying data used for the plot.

clear all

%% Step 1: Load data and region names
data_path = pwd;
load(fullfile(data_path, 'MeanStd_GraphIndex_GWM-HFN_T1.mat')) % loads results_local

% Load region information table
filename = 'E:\Neuroimage\MyProject\GMWM_Network\Data\AAL90_Descript.xlsx';
sheet = 'Name';
raw = readtable(filename, 'Sheet', sheet); % Assumes raw.ID is region index, raw.Name is label

%% Step 2: Identify hub nodes (rich-club & diversity hubs)
% Degree hubs: top 15% based on mean degree
degree_mean_matrix = results_local.DegreeSet.mean; % [90 x n_thr]
degree_column_averages = mean(degree_mean_matrix, 2); % [90 x 1]

num_elements = round(0.15 * length(degree_column_averages)); % top 15%
[~, degree_sorted_indices] = sort(degree_column_averages, 'descend');
degree_top_indices = degree_sorted_indices(1:num_elements);

% Diversity hubs: top 15% based on mean participant coefficient
participant_mean_matrix = results_local.ParticipantCoefficientSet.mean;
participant_column_averages = mean(participant_mean_matrix, 2);
[~, participant_sorted_indices] = sort(participant_column_averages, 'descend');
participant_top_indices = participant_sorted_indices(1:num_elements);

% Intersection of hubs
richclub_diversity_intersect = intersect(participant_top_indices, degree_top_indices);
zscore_hub_intersect = intersect(find(zscore(degree_column_averages)>1), degree_top_indices);

% Match region info with hub nodes and save
matched_tab = raw(ismember(raw.ID, degree_top_indices), :);
matched_tab.Degree = degree_column_averages(ismember(raw.ID, degree_top_indices));
% writetable(matched_tab, fullfile(data_path, 'HubNode.csv'));

%% Step 3: Degree distribution fitting and plotting
degree_mean_matrix = results_local.DegreeSet.mean;
degree_column_averages = mean(degree_mean_matrix, 2);

% Fit and plot degree distribution
[Para, R2] = gretna_degree_distribution(degree_column_averages, 15);

set(gca, 'FontSize', 20, 'LineWidth', 2, 'Box', 'on', 'TickDir', 'out', ...
    'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.1, 'MinorGridAlpha', 0.1, ...
    'GridLineStyle', '--', 'MinorGridLineStyle', ':', 'LineWidth', 1);


fig = gcf;
fig.Position = [1000, 600, 800, 700];

outputpath = 'E:\Neuroimage\MyProject\GMWM_Network\Code_R1\visualization\Figure4_Topo';
% print(fullfile(outputpath,"DegreeDistribution_GWM-FHN_T1.png"), '-dpng', '-r300');


%% schaefer 100
%% Step 1: Load data and region names
data_path = pwd;
load(fullfile(data_path, 'MeanStd_GraphIndex_GWM-HFN_T1_Scheafer100.mat')) % loads results_local

% Load region information table
filename = 'E:\Neuroimage\MyProject\GMWM_Network\Data\Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv';
raw = readtable(filename); % Assumes raw.ID is region index, raw.Name is label

%% Step 2: Identify hub nodes (rich-club & diversity hubs)
% Degree hubs: top 15% based on mean degree
degree_mean_matrix = results_local.DegreeSet.mean; % [90 x n_thr]
degree_column_averages = mean(degree_mean_matrix, 2); % [90 x 1]

num_elements = round(0.15 * length(degree_column_averages)); % top 15%
[~, degree_sorted_indices] = sort(degree_column_averages, 'descend');
degree_top_indices = degree_sorted_indices(1:num_elements);

% Diversity hubs: top 15% based on mean participant coefficient
participant_mean_matrix = results_local.ParticipantCoefficientSet.mean;
participant_column_averages = mean(participant_mean_matrix, 2);
[~, participant_sorted_indices] = sort(participant_column_averages, 'descend');
participant_top_indices = participant_sorted_indices(1:num_elements);

% Intersection of hubs
richclub_diversity_intersect = intersect(participant_top_indices, degree_top_indices);
zscore_hub_intersect = intersect(find(zscore(degree_column_averages)>1), degree_top_indices);

% Match region info with hub nodes and save
matched_tab = raw(ismember(raw.ROILabel, degree_top_indices), :);
matched_tab.Degree = degree_column_averages(ismember(raw.ROILabel, degree_top_indices));
% writetable(matched_tab, fullfile(data_path, 'HubNode_schaefer.csv'));

%% Step 3: Degree distribution fitting and plotting
degree_mean_matrix = results_local.DegreeSet.mean;
degree_column_averages = mean(degree_mean_matrix, 2);

% Fit and plot degree distribution
[Para, R2] = gretna_degree_distribution(degree_column_averages, 20);

set(gca, 'FontSize', 20, 'LineWidth', 2, 'Box', 'on', 'TickDir', 'out', ...
    'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.1, 'MinorGridAlpha', 0.1, ...
    'GridLineStyle', '--', 'MinorGridLineStyle', ':', 'LineWidth', 1);


fig = gcf;
fig.Position = [1000, 600, 800, 700];

outputpath = 'E:\Neuroimage\MyProject\GMWM_Network\Code_R1\visualization\Figure4_Topo';
% print(fullfile(outputpath,"DegreeDistribution_GWM-FHN_T1_schaefer.png"), '-dpng', '-r300');