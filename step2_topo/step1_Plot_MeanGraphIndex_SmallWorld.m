% Visualize small-world properties (Gamma, Lambda, Sigma) of GWM-HFN network across sparsity thresholds
% Optimized and annotated version for publication-quality figures
% Author: Zeqiang Linli, 2025-08-06

clear all

% Load the graph-theoretical attribute results for GWM-HFN (time1)
load([pwd,'\MeanStd_GraphIndex_GWM-HFN_T1.mat'])

% Define sparsity thresholds
thresholds = 0.10:0.01:0.34;

% Select the key small-world global metrics to plot
name_global = {'GammaSet', 'LambdaSet', 'SigmaSet'};

% Define colors for each metric (blue, yellow, and red)
colors = [
    0 0.4470 0.7410;      % Gamma: blue
    0.9290 0.6940 0.1250; % Lambda: yellow
    0.8500 0.3250 0.0980; % Sigma: red
];

% Custom legend labels for each curve
custom_labels = {'Gamma', 'Lambda', 'Sigma'};

%% Create figure and set properties for publication
figure;
set(gcf, 'Position', [100, 100, 800, 700]);  % Canvas size in pixels

% Preallocate handles for legend
h_lines = zeros(1, length(name_global));

% Plot each global metric (mean ± std) across thresholds
for i = 1:length(name_global)
    mean_vals = results_global.(name_global{i})(:, 1);  % Mean value across subjects
    std_vals  = results_global.(name_global{i})(:, 2);  % Standard deviation across subjects

    % Plot the mean curve
    h_lines(i) = plot(thresholds, mean_vals, 'LineWidth', 2, 'Color', colors(i, :)); hold on;

    % Add shaded area for ± 1 SD (confidence band)
    fill([thresholds, fliplr(thresholds)], ...
         [mean_vals - std_vals; flipud(mean_vals + std_vals)], ...
         colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% Add reference horizontal line at y = 1 (small-world criterion)
yline(1, '--', 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840]);

% Add legend for mean curves (ignore fill bands)
legend(h_lines, custom_labels, 'FontSize', 15, 'FontWeight','bold','Location', 'northeast');

% Axis labels
xlabel('Sparsity', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Value',    'FontSize', 20, 'FontWeight', 'bold');

% Axis limits for clarity
xlim([0.10, 0.35]);
ylim([0, 7]);

% Optionally auto-adjust y-limits based on data (uncomment if desired)
% y_min = min(cellfun(@(x) min(results_global.(x)(:, 1) - results_global.(x)(:, 2)), name_global));
% y_max = max(cellfun(@(x) max(results_global.(x)(:, 1) + results_global.(x)(:, 2)), name_global));
% ylim([y_min, y_max]);

% Beautify plot: grid, font, box, etc.
set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'Box', 'on', 'TickDir', 'out', ...
    'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.1, 'MinorGridAlpha', 0.1, ...
    'GridLineStyle', '--', 'MinorGridLineStyle', ':', 'LineWidth', 1);

% Save figure as high-resolution TIFF for publication
% print(gcf, 'MeanGraphIndex_Sigma_Lambda_Gamma_GWM-FHN_T1', '-dtiffn', '-r300');

% End of script



%% schaefer 100

% Load the graph-theoretical attribute results for GWM-HFN (time1)
load('E:\Neuroimage\MyProject\GMWM_Network\Code_R1\step2_topo\MeanStd_GraphIndex_GWM-HFN_T1_Scheafer100.mat')

% Define sparsity thresholds
thresholds = 0.10:0.01:0.34;

% Select the key small-world global metrics to plot
name_global = {'GammaSet', 'LambdaSet', 'SigmaSet'};

% Define colors for each metric (blue, yellow, and red)
colors = [
    0 0.4470 0.7410;      % Gamma: blue
    0.9290 0.6940 0.1250; % Lambda: yellow
    0.8500 0.3250 0.0980; % Sigma: red
];

% Custom legend labels for each curve
custom_labels = {'Gamma', 'Lambda', 'Sigma'};

%% Create figure and set properties for publication
figure;
set(gcf, 'Position', [100, 100, 800, 700]);  % Canvas size in pixels

% Preallocate handles for legend
h_lines = zeros(1, length(name_global));

% Plot each global metric (mean ± std) across thresholds
for i = 1:length(name_global)
    mean_vals = results_global.(name_global{i})(:, 1);  % Mean value across subjects
    std_vals  = results_global.(name_global{i})(:, 2);  % Standard deviation across subjects

    % Plot the mean curve
    h_lines(i) = plot(thresholds, mean_vals, 'LineWidth', 2, 'Color', colors(i, :)); hold on;

    % Add shaded area for ± 1 SD (confidence band)
    fill([thresholds, fliplr(thresholds)], ...
         [mean_vals - std_vals; flipud(mean_vals + std_vals)], ...
         colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% Add reference horizontal line at y = 1 (small-world criterion)
yline(1, '--', 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840]);

% Add legend for mean curves (ignore fill bands)
legend(h_lines, custom_labels, 'FontSize', 15, 'FontWeight','bold','Location', 'northeast');

% Axis labels
xlabel('Sparsity', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Value',    'FontSize', 20, 'FontWeight', 'bold');

% Axis limits for clarity
xlim([0.10, 0.35]);
ylim([0, 7]);

% Optionally auto-adjust y-limits based on data (uncomment if desired)
% y_min = min(cellfun(@(x) min(results_global.(x)(:, 1) - results_global.(x)(:, 2)), name_global));
% y_max = max(cellfun(@(x) max(results_global.(x)(:, 1) + results_global.(x)(:, 2)), name_global));
% ylim([y_min, y_max]);

% Beautify plot: grid, font, box, etc.
set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'Box', 'on', 'TickDir', 'out', ...
    'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.1, 'MinorGridAlpha', 0.1, ...
    'GridLineStyle', '--', 'MinorGridLineStyle', ':', 'LineWidth', 1);

% Save figure as high-resolution TIFF for publication
% print(gcf, 'MeanGraphIndex_Sigma_Lambda_Gamma_GWM-FHN_T1_schaefer', '-dtiffn', '-r300');

% End of script