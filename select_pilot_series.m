function pilot_series = select_pilot_series(segments_pre, segments_post, intensity_levels)
%% Randomly select one preTMS and one postTMS segment from each intensity level
% Inputs:
%   segments_pre - structure with preTMS segments
%   segments_post - structure with postTMS segments
%   intensity_levels - cell array of intensity level names
% Outputs:
%   pilot_series - cell array of 6 time series (3 preTMS + 3 postTMS)

    pilot_series = cell(6, 1);
    series_labels = cell(6, 1);
    
    for i = 1:length(intensity_levels)
        level = intensity_levels{i};
        
        % Get preTMS data (standardized: [channels x time x episodes])
        pre_data = segments_pre.(level);
        dims_pre = size(pre_data);
        if length(dims_pre) ~= 3
            error('Expected 3D data for preTMS %s', level);
        end
        n_channels     = dims_pre(1);
        n_episodes_pre = dims_pre(3);
        
        % Get postTMS data (standardized: [channels x time x episodes])
        post_data = segments_post.(level);
        dims_post = size(post_data);
        if length(dims_post) ~= 3
            error('Expected 3D data for postTMS %s', level);
        end
        n_episodes_post = dims_post(3);
        
        % Randomly select one episode from each
        selected_episode_pre = randi(n_episodes_pre);
        selected_episode_post = randi(n_episodes_post);
        
        % Extract time series (use first channel if multiple channels)
        pilot_series{i}   = squeeze(pre_data(1, :, selected_episode_pre));
        pilot_series{i+3} = squeeze(post_data(1, :, selected_episode_post));
        
        series_labels{i} = sprintf('%s_preTMS_ep%d', level, selected_episode_pre);
        series_labels{i+3} = sprintf('%s_postTMS_ep%d', level, selected_episode_post);
        
        fprintf('Selected pilot series %d: %s (length: %d)\n', ...
            i, series_labels{i}, length(pilot_series{i}));
        fprintf('Selected pilot series %d: %s (length: %d)\n', ...
            i+3, series_labels{i+3}, length(pilot_series{i+3}));
    end
    
    % Store labels in a persistent way (could use a structure)
    assignin('caller', 'pilot_series_labels', series_labels);
end
