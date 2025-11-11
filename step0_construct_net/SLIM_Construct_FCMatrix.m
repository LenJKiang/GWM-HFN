% This script computes and saves GM-GM and GWM-HFN matrices
% for both all subjects (per timepoint) and for matched subjects (common IDs) 
% in both time1 and time2 of the SLIM dataset.
% Author: Linli Ze-Qiang, 2025-07-23


clear; clc;
%% Path settings 
dir_time1 = 'E:\Neuroimage\MyProject\GMWM_Network\Data\SLIM\time1_ROISignal';
dir_time2 = 'E:\Neuroimage\MyProject\GMWM_Network\Data\SLIM\time2_ROISignal';
save_dir = 'E:\Neuroimage\MyProject\GMWM_Network\Data\SLIM\'; 

% List ROI signal files
files_time1 = dir(fullfile(dir_time1, 'ROISignals_sub-*.txt'));
files_time2 = dir(fullfile(dir_time2, 'ROISignals_sub-*.txt'));
% Extract subject IDs (strings)
sub_time1 = extractBetween({files_time1.name}, 'ROISignals_', '.txt');
sub_time2 = extractBetween({files_time2.name}, 'ROISignals_', '.txt');

%% Parameters
nGM = 90; nWM = 48;

%% 1. Compute for all time1 subjects
n1 = numel(files_time1);
GMGM_time1_all = zeros(nGM, nGM, n1);
GWM_HFN_time1_all = zeros(nGM, nGM, n1);

for i = 1:n1
    ts = load(fullfile(dir_time1, files_time1(i).name)); % [T x 138]
    GM = ts(:,1:nGM);
    WM = ts(:,nGM+1:end);
    GMGM_time1_all(:,:,i) = corr(GM);
    GMWM = corr(GM, WM);
    zGMWM = zscore(GMWM,0,2); % row-wise
    GWM_HFN_time1_all(:,:,i) = zGMWM * zGMWM';
end

%% 2. Compute for all time2 subjects
n2 = numel(files_time2);
GMGM_time2_all = zeros(nGM, nGM, n2);
GWM_HFN_time2_all = zeros(nGM, nGM, n2);

for i = 1:n2
    ts = load(fullfile(dir_time2, files_time2(i).name)); % [T x 138]
    GM = ts(:,1:nGM);
    WM = ts(:,nGM+1:end);
    GMGM_time2_all(:,:,i) = corr(GM);
    GMWM = corr(GM, WM);
    zGMWM = zscore(GMWM,0,2);
    GWM_HFN_time2_all(:,:,i) = zGMWM * zGMWM';
end

%% 3. Find common subjects (intersection, keep order in time1 and time2)
[common_subs, idx_time1, idx_time2] = intersect(sub_time1, sub_time2, 'stable');
n_common = numel(common_subs);

% Create "common" subsets (retain same naming pattern: time{t}_<NET>_common)
time1_GMGM_common = GMGM_time1_all(:,:,idx_time1);
time1_GWM_HFN_common = GWM_HFN_time1_all(:,:,idx_time1);
time2_GMGM_common = GMGM_time2_all(:,:,idx_time2);
time2_GWM_HFN_common = GWM_HFN_time2_all(:,:,idx_time2);

% Also produce variables following unified pattern for "all" datasets:
time1_GMGM = GMGM_time1_all;
time1_GWM_HFN = GWM_HFN_time1_all;
time2_GMGM = GMGM_time2_all;
time2_GWM_HFN = GWM_HFN_time2_all;

%% 4. Save all results (with subject IDs) using unified variable names
cd(save_dir);

% time{t}_NET.mat for ALL subjects (variable: time{t}_NET)
save('GMGM_time1.mat','time1_GMGM','sub_time1','-v6');
save('GMGM_time2.mat','time2_GMGM','sub_time2','-v6');
save('GWM_HFN_time1.mat','time1_GWM_HFN','sub_time1','-v6');
save('GWM_HFN_time2.mat','time2_GWM_HFN','sub_time2','-v6');

% time{t}_NET_common.mat for COMMON matched subjects (variable: time{t}_NET_common)
save('GMGM_time1_common.mat','time1_GMGM_common','common_subs','-v6');
save('GMGM_time2_common.mat','time2_GMGM_common','common_subs','-v6');
save('GWM_HFN_time1_common.mat','time1_GWM_HFN_common','common_subs','-v6');
save('GWM_HFN_time2_common.mat','time2_GWM_HFN_common','common_subs','-v6');

% Also save subject id lists in loader-friendly names time{t}_subject_ids.mat (variable: time{t}_subject_ids)
time1_subject_ids = sub_time1;
time2_subject_ids = sub_time2;
save('time1_subject_ids.mat','time1_subject_ids','-v6');
save('time2_subject_ids.mat','time2_subject_ids','-v6');

disp('All network matrices saved for both all and common matched subjects with unified variable names.');