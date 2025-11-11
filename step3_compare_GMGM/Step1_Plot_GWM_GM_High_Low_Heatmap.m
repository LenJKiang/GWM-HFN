% Visualization of Low and High Correlation Edges between GWM-HFN and GM-GM
% Optimization: English comments, enhanced colorbar font size and label position

output_path = 'E:\Neuroimage\MyProject\GMWM_Network\Code_R1\visualization\Figure5_CompareGMGM'

load("GWM_GM_HighCorrEdges.mat")
load("GWM_GM_LowCorrEdges.mat")

% Load network indices and legend labels
netIndex = readmatrix("E:\Neuroimage\MyProject\GMWM_Network\Data\AAL_7networkIndex.csv");
legends = {'VN', 'SMN', 'AN', 'LB', 'FPN', 'DMN', 'BGN'};

%% Visualize Low-Correlation Edges
figure('Position', [100, 100, 300, 300]); % Keep original figure size

LowCorrMatrix(LowCorrMatrix == 0) = 1; % Set zeros to ones for background

lc_netplot('-n', LowCorrMatrix, '-ni', netIndex, '-il', 1, '-lg', legends);
axis square  % Make the plot square for better aesthetics

colormap("hot");

colorbar_handle = colorbar;
colorbar_handle.Location = 'north'; % Place colorbar at the top
colorbar_handle.Position = [0.19, 0.85, 0.65, 0.03]; % Custom position
colorbar_handle.TickDirection = 'in'; % Ticks inward

caxis([0.60, 0.75]);  % Color axis limits
colorbar_handle.Ticks = 0.6:0.05:0.75;

colorbar_handle.FontSize = 10; % Increase font size for publication

% saveas(gcf, fullfile(output_path,'GWM_GM_LowCorrEdges.pdf'));  % Save as PDF for high-quality output

%% Visualize High-Correlation Edges
figure('Position', [100, 100, 300, 300]);

HighCorrMatrix(HighCorrMatrix == 0) = 1;

lc_netplot('-n', HighCorrMatrix, '-ni', netIndex, '-il', 1, '-lg', legends);
axis square

colormap("hot");

colorbar_handle = colorbar;
colorbar_handle.Location = 'north';
colorbar_handle.Position = [0.19, 0.85, 0.65, 0.03];
colorbar_handle.TickDirection = 'in';

caxis([0.80, 0.95]);
colorbar_handle.Ticks = 0.8:0.05:0.95;

colorbar_handle.FontSize = 10; % Increase font size

% saveas(gcf, fullfile(output_path,'GWM_GM_HighCorrEdges.pdf'));  % Save as PDF for high-quality output

