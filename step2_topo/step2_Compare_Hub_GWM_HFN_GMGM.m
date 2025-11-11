% This script compares hub regions between GWM-HFN (GM-WM-GM) and traditional GM-GM networks.
% 1. Extracts hub nodes (top 15% by degree) from both methods.
% 2. Performs quantitative overlap analysis (Dice coefficient, Jaccard index).
% Author: Zeqiang LINLI, 2025-08-06

clear all; clc;

%% Step 1: Load degree data for both networks
% --- GWM-HFN (GM-WM-GM) ---
GMWMGM = load([pwd,'\MeanStd_GraphIndex_GWM-HFN_T1.mat']);
degree_gmwmgm = mean(GMWMGM.results_local.DegreeSet.mean, 2); % (90x1)

% --- GM-GM ---
GMGM = load([pwd,'\MeanStd_GraphIndex_GMGM_T1.mat']);
degree_gmgm = mean(GMGM.results_local.DegreeSet.mean, 2); % (90x1)

%% Step 2: Extract hub nodes (top 15%)
num_nodes = length(degree_gmwmgm);
num_hub = round(0.15 * num_nodes);

% GWM-HFN hubs
[~, idx_sorted_gmwmgm] = sort(degree_gmwmgm, 'descend');
hub_gmwmgm = idx_sorted_gmwmgm(1:num_hub);

% GM-GM hubs
[~, idx_sorted_gmgm] = sort(degree_gmgm, 'descend');
hub_gmgm = idx_sorted_gmgm(1:num_hub);

%% Step 3: Quantitative overlap analysis
% Dice coefficient: 2 * |A ∩ B| / (|A| + |B|)
n_common = numel(intersect(hub_gmwmgm, hub_gmgm));
dice_coef = 2 * n_common / (num_hub + num_hub);

% Jaccard index: |A ∩ B| / |A ∪ B|
jaccard_index = n_common / numel(union(hub_gmwmgm, hub_gmgm));

fprintf('Quantitative comparison of hub regions:\n');
fprintf('  Number of hubs (each method): %d\n', num_hub);
fprintf('  Common hubs: %d\n', n_common);
fprintf('  Dice coefficient: %.3f\n', dice_coef);
fprintf('  Jaccard index: %.3f\n', jaccard_index);

%% Step 4: Save results for reporting
results_hub_overlap = struct();
results_hub_overlap.hub_gmwmgm = hub_gmwmgm;
results_hub_overlap.hub_gmgm = hub_gmgm;
results_hub_overlap.common_hub = intersect(hub_gmwmgm, hub_gmgm);
results_hub_overlap.gmwmgm_specific = setdiff(hub_gmwmgm, hub_gmgm);
results_hub_overlap.gmgm_specific = setdiff(hub_gmgm, hub_gmwmgm);
results_hub_overlap.dice_coef = dice_coef;
results_hub_overlap.jaccard_index = jaccard_index;

% save('results_hub_overlap.mat', 'results_hub_overlap');

