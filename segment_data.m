function [segments_pre, segments_post] = segment_data(data_clean, sampling_freq, pre_duration_ms, post_duration_ms)
%% Segment data into preTMS (first 900ms) and postTMS (last 900ms) periods
% Inputs:
%   data_clean - structure with clean episodes for each intensity level
%   sampling_freq - sampling frequency in Hz
%   pre_duration_ms - duration of preTMS period in ms
%   post_duration_ms - duration of postTMS period in ms
% Outputs:
%   segments_pre - structure with preTMS segments for each intensity level
%   segments_post - structure with postTMS segments for each intensity level

    segments_pre = struct();
    segments_post = struct();
    intensity_levels = fieldnames(data_clean);
    
    % Convert ms to samples
    pre_samples = round(pre_duration_ms / 1000 * sampling_freq);
    post_samples = round(post_duration_ms / 1000 * sampling_freq);
    
    fprintf('Segmenting: preTMS = %d samples (%.1f ms), postTMS = %d samples (%.1f ms)\n', ...
        pre_samples, pre_duration_ms, post_samples, post_duration_ms);
    
    for i = 1:length(intensity_levels)
        level = intensity_levels{i};
        data = data_clean.(level);
        
        % Determine data dimensions
        dims = size(data);
        if length(dims) == 3
            if dims(2) > dims(1) && dims(2) > dims(3)
                % Format: [channels x time x episodes]
                n_channels = dims(1);
                n_samples = dims(2);
                n_episodes = dims(3);
                time_dim = 2;
            elseif dims(3) > dims(1) && dims(3) > dims(2)
                % Format: [episodes x channels x time]
                n_episodes = dims(1);
                n_channels = dims(2);
                n_samples = dims(3);
                time_dim = 3;
            else
                error('Cannot determine data dimensions for %s', level);
            end
        else
            error('Expected 3D data for %s', level);
        end
        
        if n_samples < pre_samples + post_samples
            warning('Insufficient samples for %s: have %d, need at least %d', ...
                level, n_samples, pre_samples + post_samples);
        end
        
        % Extract preTMS (first pre_samples) and postTMS (last post_samples)
        if time_dim == 2
            % Format: [channels x time x episodes]
            segments_pre.(level) = data(:, 1:pre_samples, :);
            segments_post.(level) = data(:, (end-post_samples+1):end, :);
        else
            % Format: [episodes x channels x time]
            segments_pre.(level) = data(:, :, 1:pre_samples);
            segments_post.(level) = data(:, :, (end-post_samples+1):end);
        end
        
        fprintf('Segmented %s: preTMS size %s, postTMS size %s\n', ...
            level, mat2str(size(segments_pre.(level))), ...
            mat2str(size(segments_post.(level))));
    end
end
