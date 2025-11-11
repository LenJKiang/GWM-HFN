% This script constructs GM-GM, and GWM_HFN connectivity matrices
% for all subjects at three sessions in the BNU-3 dataset. 
% Output variables follow the unified pattern: time{t}_<NET> (e.g. time1_GMGM)
% Author: Linli Ze-Qiang, adapted 2025-07-23

clear; clc;

% Define data root and session folders
data_root = 'E:\Neuroimage\MyProject\GMWM_Network\Data\BNU3\'; 
session_folders = {'S1_ROISignals', 'S2_ROISignals', 'S3_ROISignals'};

% Loop over sessions
for s = 1:3
    fprintf('\n==== Processing Session %d ====\n', s);
    session_dir = fullfile(data_root, session_folders{s});
    files = dir(fullfile(session_dir, 'ROISignals_sub-*.txt'));
    num_subjects = numel(files);
    fprintf('Found %d subjects in %s\n', num_subjects, session_dir);

    % Preallocate matrices
    time_var_GMGM = zeros(90, 90, num_subjects);
    time_var_WMWM = zeros(48, 48, num_subjects);
    time_var_GMWM = zeros(90, 48, num_subjects);
    time_var_GWM_HFN = zeros(90, 90, num_subjects);

    % Loop over subjects
    for i = 1:num_subjects
        % Load ROI signals (T x 138)
        data = load(fullfile(session_dir, files(i).name));
        if size(data,2) ~= 138
            error('File %s: ROI signal columns mismatch (should be 138)', files(i).name);
        end
        sig_gm = data(:,1:90);   % Gray matter (GM)
        sig_wm = data(:,91:138); % White matter (WM)

        % GM-GM functional connectivity
        time_var_GMGM(:, :, i) = corrcoef(sig_gm);
        % WM-WM functional connectivity
        time_var_WMWM(:, :, i) = corrcoef(sig_wm);
        % GM-WM functional connectivity
        time_var_GMWM(:, :, i) = corr(sig_gm, sig_wm);
    end

    % Compute GWM-HFN matrix for each subject (row-normalized GMWM times itself)
    for i = 1:num_subjects
        GMWM = time_var_GMWM(:,:,i);
        zGMWM = zscore(GMWM'); % 48 x 90, row-normalize
        time_var_GWM_HFN(:,:,i) = zGMWM' * zGMWM; % 90x90
    end

    % Assign to unified variable names: time{s}_NET
    switch s
        case 1
            time1_GMGM = time_var_GMGM;
            time1_WMWM = time_var_WMWM;
            time1_GMWM = time_var_GMWM;
            time1_GWM_HFN = time_var_GWM_HFN;
        case 2
            time2_GMGM = time_var_GMGM;
            time2_WMWM = time_var_WMWM;
            time2_GMWM = time_var_GMWM;
            time2_GWM_HFN = time_var_GWM_HFN;
        case 3
            time3_GMGM = time_var_GMGM;
            time3_WMWM = time_var_WMWM;
            time3_GMWM = time_var_GMWM;
            time3_GWM_HFN = time_var_GWM_HFN;
    end

    % Save results as files named <NET>_time{t}.mat with variable time{t}_<NET>
    % e.g. GMGM_time1.mat contains variable time1_GMGM
    save(fullfile(data_root, sprintf('GMGM_time%d.mat', s)), sprintf('time%d_GMGM', s), '-v6');
    save(fullfile(data_root, sprintf('GWM_HFN_time%d.mat', s)), sprintf('time%d_GWM_HFN', s), '-v6');
    save(fullfile(data_root, sprintf('GMWM_time%d.mat', s)), sprintf('time%d_GMWM', s), '-v6');
    save(fullfile(data_root, sprintf('WMWM_time%d.mat', s)), sprintf('time%d_WMWM', s), '-v6');

    % Save subject ids for this session (use time{t}_subject_ids variable)
    ids = cell(num_subjects,1);
    for i = 1:num_subjects
        tokens = regexp(files(i).name, 'ROISignals_sub-(\w+).txt', 'tokens');
        if ~isempty(tokens)
            ids{i} = tokens{1}{1};
        else
            ids{i} = files(i).name;
        end
    end
    eval(sprintf('time%d_subject_ids = ids;', s));
    save(fullfile(data_root, sprintf('time%d_subject_ids.mat', s)), sprintf('time%d_subject_ids', s), '-v6');

    fprintf('Session %d done: matrices saved with unified variable names.\n', s);
end