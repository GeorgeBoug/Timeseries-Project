%% Main Workflow Script for TMS-EEG Time Series Analysis
% This script orchestrates the complete analysis pipeline:
% 1. Load EEG data and channel names
% 2. Select electrodes based on AEM numbers
% 3. Visualize episodes for manual rejection
% 4. Segment data into preTMS and postTMS periods
% 5. Pilot analysis on 6 series (linear and nonlinear measures)
% 6. Full analysis on all segments
% 7. Statistical comparisons and visualization

clear; close all; clc;

%% Add lab function folders to path (avoid name conflicts)
% NOTE: both lab folders include some functions with the same name
% (e.g. 'autocorrelation', 'nrmse'). We add labnonlinear first and
% lablinear second so that lablinear versions take precedence when names
% collide.
thisFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(thisFolder, 'labnonlinear'));
addpath(fullfile(thisFolder, 'lablinear'));

%% Configuration
AEM_NUMBERS = [1, 1];  % Team members' AEM numbers
SAMPLING_FREQ = 1450;  % Hz
PRE_TMS_DURATION = 900;  % ms
POST_TMS_DURATION = 900;  % ms
N_CLEAN_EPISODES = 20;  % Number of clean episodes per intensity level
PILOT_SAMPLES_FIT = 600;  % Samples for model fitting
PILOT_SAMPLES_FORECAST = 300;  % Samples for forecasting
FORECAST_HORIZONS = 1:10;  % Step-ahead forecast horizons

% File paths
DATA_FOLDER = '.';  % Current folder (change to 'data' if files are in subfolder)
CHANNELS_FILE = 'Channels.txt';

% Auto-detect .mat files and map to intensity levels
fprintf('Auto-detecting .mat files...\n');
mat_files = dir(fullfile(DATA_FOLDER, '*.mat'));
if isempty(mat_files)
    error('No .mat files found in %s', DATA_FOLDER);
end

% Extract intensity values from filenames (e.g., I036 -> 36%, I054 -> 54%, I072 -> 72%)
intensity_values = zeros(length(mat_files), 1);
for i = 1:length(mat_files)
    filename = mat_files(i).name;
    % Find pattern I### in filename
    match = regexp(filename, 'I(\d+)', 'tokens');
    if ~isempty(match)
        intensity_values(i) = str2double(match{1}{1});
    else
        intensity_values(i) = i * 30;  % Fallback: assign sequential values
    end
end

% Sort by intensity and assign labels
[~, sort_idx] = sort(intensity_values);
mat_files = mat_files(sort_idx);
INTENSITY_LEVELS = {'low', 'medium', 'high'};

if length(mat_files) < 3
    error('Found only %d .mat files, need 3 files for three intensity levels', length(mat_files));
end

% Use only the first 3 files (sorted by intensity)
MAT_FILE_PATHS = cell(3, 1);
for i = 1:3
    MAT_FILE_PATHS{i} = fullfile(DATA_FOLDER, mat_files(i).name);
    fprintf('  Intensity %s: %s (I%d)\n', INTENSITY_LEVELS{i}, mat_files(i).name, intensity_values(sort_idx(i)));
end

%% Step 1: Load data
fprintf('\nStep 1: Loading data...\n');
[data_all, channel_names] = load_eeg_data(MAT_FILE_PATHS, CHANNELS_FILE, INTENSITY_LEVELS);

%% Step 2: Select electrodes
fprintf('Step 2: Selecting electrodes...\n');
[electrode_indices, electrode_names] = select_electrodes(AEM_NUMBERS, channel_names);
fprintf('Selected electrodes: %s (index %d) and %s (index %d)\n', ...
    electrode_names{1}, electrode_indices(1), electrode_names{2}, electrode_indices(2));

%% Step 3: Extract selected channels
fprintf('Step 3: Extracting selected channels...\n');
data_selected = extract_channels(data_all, electrode_indices);

%% Step 4: Visualize episodes for manual rejection
fprintf('Step 4: Visualizing episodes for manual rejection...\n');
rejected_episodes = visualize_and_reject_episodes(data_selected, INTENSITY_LEVELS, SAMPLING_FREQ);

%% Step 5: Select clean episodes
fprintf('Step 5: Selecting clean episodes...\n');
data_clean = select_clean_episodes(data_selected, rejected_episodes, N_CLEAN_EPISODES);

%% Step 6: Segment data into preTMS and postTMS
fprintf('Step 6: Segmenting data...\n');
[segments_pre, segments_post] = segment_data(data_clean, SAMPLING_FREQ, ...
    PRE_TMS_DURATION, POST_TMS_DURATION);

%% Step 7: Select pilot series (one preTMS and one postTMS from each intensity level)
fprintf('Step 7: Selecting pilot series...\n');
pilot_series = select_pilot_series(segments_pre, segments_post, INTENSITY_LEVELS);

%% Step 8: Pilot linear analysis
fprintf('Step 8: Pilot linear analysis...\n');
[acf_pilot, pacf_pilot, nrmse_pilot, linear_measure_name, linear_measure_params] = ...
    pilot_linear_analysis(pilot_series, SAMPLING_FREQ, PILOT_SAMPLES_FIT, ...
    PILOT_SAMPLES_FORECAST, FORECAST_HORIZONS);

%% Step 9: Pilot nonlinear analysis
fprintf('Step 9: Pilot nonlinear analysis...\n');
[nonlinear_measures_pilot, nonlinear_measure_name, nonlinear_measure_params] = ...
    pilot_nonlinear_analysis(pilot_series, SAMPLING_FREQ);

%% Step 10: Full analysis - compute measures for all segments
fprintf('Step 10: Full analysis on all segments...\n');
[linear_values_pre, linear_values_post] = compute_linear_measures_all(...
    segments_pre, segments_post, linear_measure_name, linear_measure_params, ...
    SAMPLING_FREQ, PILOT_SAMPLES_FIT, PILOT_SAMPLES_FORECAST, FORECAST_HORIZONS);

[nonlinear_values_pre, nonlinear_values_post] = compute_nonlinear_measures_all(...
    segments_pre, segments_post, nonlinear_measure_name, nonlinear_measure_params, SAMPLING_FREQ);

%% Step 11: Statistical comparisons and visualization
fprintf('Step 11: Statistical comparisons and visualization...\n');
perform_statistical_analysis(linear_values_pre, linear_values_post, ...
    nonlinear_values_pre, nonlinear_values_post, INTENSITY_LEVELS, ...
    linear_measure_name, nonlinear_measure_name);

fprintf('Analysis complete!\n');
