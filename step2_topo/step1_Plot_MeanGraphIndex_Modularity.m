% Visualize modularity (Q value) of the GWM-HFN network across sparsity thresholds
% Optimized and annotated for publication-quality figure
% Author: Zeqiang Linli, 2025-08-06

clear all

% Set working directory for saving figures
% Load the modularity results for GWM-HFN (time1)
load([pwd,'\MeanStd_GraphIndex_GWM-HFN_T1.mat'])

% Define sparsity thresholds
thresholds = 0.10:0.01:0.34;

% Select the modularity global metric to plot
name_global = {'ModularitySet'};

% Define color for modularity curve (red)
colors = [
    0.8500 0.3250 0.0980;  % Modularity: red
];

% Custom legend label
custom_labels = {'Modularity'};

%% Create figure 
% AAL90
figure;
set(gcf, 'Position', [100, 100, 800, 700]);  % Canvas size in pixels

% Preallocate handle for legend
h_lines = zeros(1, length(name_global));

% Plot modularity (mean ± std) across thresholds
for i = 1:length(name_global)
    mean_vals = results_global.(name_global{i})(:, 1);  % Mean modularity Q value
    std_vals  = results_global.(name_global{i})(:, 2);  % Standard deviation

    % Plot the mean curve
    h_lines(i) = plot(thresholds, mean_vals, 'LineWidth', 2, 'Color', colors(i, :)); hold on;

    % Add shaded area for ± 1 SD (confidence band)
    fill([thresholds, fliplr(thresholds)], ...
         [mean_vals - std_vals; flipud(mean_vals + std_vals)], ...
         colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% Add reference horizontal line at Q = 0.30 (modularity threshold)
yline(0.30, '--', 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840]);

% Add legend for mean curve (ignore fill band)
legend(h_lines, custom_labels, 'FontSize', 15, 'FontWeight','bold','Location', 'northeast');

% Axis labels
xlabel('Sparsity', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Q Value', 'FontSize', 20, 'FontWeight', 'bold');

% Axis limits for clarity
xlim([0.10, 0.35]);
ylim([0.20, 0.6]);

% Beautify plot: grid, font, box, etc.
set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'Box', 'on', 'TickDir', 'out', ...
    'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.1, 'MinorGridAlpha', 0.1, ...
    'GridLineStyle', '--', 'MinorGridLineStyle', ':', 'LineWidth', 1);

% Save figure as high-resolution TIFF for publication
% print(gcf, 'MeanGraphIndex_Modularity_GWM-FHN_T1', '-dtiffn', '-r300');

% End of script



%% schaefer 100
% Load the modularity results for GWM-HFN (time1)
load([pwd,'\MeanStd_GraphIndex_GWM-HFN_T1_Scheafer100.mat'])
% Define sparsity thresholds
thresholds = 0.10:0.01:0.34;

% Select the modularity global metric to plot
name_global = {'ModularitySet'};

% Define color for modularity curve (red)
colors = [
    0.8500 0.3250 0.0980;  % Modularity: red
];

% Custom legend label
custom_labels = {'Modularity'};

%% Create figure and set properties for publication
figure;
set(gcf, 'Position', [100, 100, 800, 700]);  % Canvas size in pixels

% Preallocate handle for legend
h_lines = zeros(1, length(name_global));

% Plot modularity (mean ± std) across thresholds
for i = 1:length(name_global)
    mean_vals = results_global.(name_global{i})(:, 1);  % Mean modularity Q value
    std_vals  = results_global.(name_global{i})(:, 2);  % Standard deviation

    % Plot the mean curve
    h_lines(i) = plot(thresholds, mean_vals, 'LineWidth', 2, 'Color', colors(i, :)); hold on;

    % Add shaded area for ± 1 SD (confidence band)
    fill([thresholds, fliplr(thresholds)], ...
         [mean_vals - std_vals; flipud(mean_vals + std_vals)], ...
         colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% Add reference horizontal line at Q = 0.30 (modularity threshold)
yline(0.30, '--', 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840]);

% Add legend for mean curve (ignore fill band)
legend(h_lines, custom_labels, 'FontSize', 15, 'FontWeight','bold','Location', 'northeast');

% Axis labels
xlabel('Sparsity', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Q Value', 'FontSize', 20, 'FontWeight', 'bold');

% Axis limits for clarity
xlim([0.10, 0.35]);
ylim([0.20, 0.6]);

% Beautify plot: grid, font, box, etc.
set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'Box', 'on', 'TickDir', 'out', ...
    'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.1, 'MinorGridAlpha', 0.1, ...
    'GridLineStyle', '--', 'MinorGridLineStyle', ':', 'LineWidth', 1);

% Save figure as high-resolution TIFF for publication
% print(gcf, 'MeanGraphIndex_Modularity_GWM-FHN_T1_Schaefer', '-dtiffn', '-r300');

% End of script


