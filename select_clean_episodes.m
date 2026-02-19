function data_clean = select_clean_episodes(data_selected, rejected_episodes, n_clean)
%% Select N clean episodes from each intensity level
% Inputs:
%   data_selected - structure with selected channels for each intensity level
%   rejected_episodes - structure with rejected episode indices
%   n_clean - number of clean episodes to select per intensity level
% Outputs:
%   data_clean - structure with clean episodes for each intensity level

    data_clean = struct();
    intensity_levels = fieldnames(data_selected);
    
    for i = 1:length(intensity_levels)
        level = intensity_levels{i};
        data = data_selected.(level);
        rejected = rejected_episodes.(level);
        
        % After extraction, data is standardized to [channels x time x episodes]
        dims = size(data);
        if length(dims) ~= 3
            error('Expected 3D data for %s', level);
        end
        n_episodes   = dims(3);
        
        % Get clean episode indices
        all_episodes = 1:n_episodes;
        clean_episodes = setdiff(all_episodes, rejected);
        
        if length(clean_episodes) < n_clean
            warning('Only %d clean episodes available for %s, requested %d. Using all available.', ...
                length(clean_episodes), level, n_clean);
            selected_indices = clean_episodes;
        else
            % Randomly select n_clean episodes
            selected_indices = clean_episodes(randperm(length(clean_episodes), n_clean));
            selected_indices = sort(selected_indices);
        end
        
        % Extract selected episodes (keep [channels x time x episodes])
        data_clean.(level) = data(:, :, selected_indices);
        
        fprintf('Selected %d clean episodes for %s (indices: %s)\n', ...
            length(selected_indices), level, mat2str(selected_indices));
    end
end
