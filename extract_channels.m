function data_selected = extract_channels(data_all, electrode_indices)
%% Extract only the selected channels from all intensity level datasets
% Inputs:
%   data_all - structure with cMF arrays for each intensity level
%   electrode_indices - array of electrode indices to extract
% Outputs:
%   data_selected - structure with selected channels for each intensity level

    data_selected = struct();
    intensity_levels = fieldnames(data_all);
    
    for i = 1:length(intensity_levels)
        level = intensity_levels{i};
        cMF_data = data_all.(level);
        
        % According to the dataset description, cMF has dimensions:
        %   n x K x M  (time x channels x episodes), where
        %   n = 2901 samples, K = 60 channels, M = number of episodes.
        %
        % Here we:
        %   1) Detect which dimension is channels and which is time
        %   2) Extract the requested electrodes along the channel dimension
        %   3) Standardize the result to [channels x time x episodes]
        
        dims = size(cMF_data);
        if length(dims) ~= 3
            error('Unexpected data dimensions for %s: %s', level, mat2str(dims));
        end
        
        % --- Detect channel dimension ---
        % Prefer a dimension of size 60; fall back to any dimension large
        % enough to index all requested electrodes.
        chan_dim = find(dims == 60, 1);
        if isempty(chan_dim)
            candidates = find(dims >= max(electrode_indices));
            if isempty(candidates)
                error('Cannot determine channel dimension for %s', level);
            end
            % Choose the smallest suitable dimension as channels
            [~, idx_min] = min(dims(candidates));
            chan_dim = candidates(idx_min);
        end
        
        % --- Detect time dimension ---
        % Prefer the dimension with 2901 samples; otherwise take the largest
        % remaining dimension.
        time_dim = find(dims == 2901, 1);
        if isempty(time_dim)
            other_dims = setdiff(1:3, chan_dim);
            [~, idx_max] = max(dims(other_dims));
            time_dim = other_dims(idx_max);
        end
        
        % The remaining dimension is episodes
        ep_dim = setdiff(1:3, [chan_dim, time_dim]);
        if numel(ep_dim) ~= 1
            error('Cannot determine episode dimension for %s', level);
        end
        ep_dim = ep_dim(1);
        
        % --- Extract requested electrodes along the channel dimension ---
        idx = {':', ':', ':'};
        idx{chan_dim} = electrode_indices;
        sub_data = cMF_data(idx{:});
        
        % --- Standardize to [channels x time x episodes] ---
        % Current order is [d1 d2 d3] with known roles (chan_dim, time_dim, ep_dim)
        perm_order = [chan_dim, time_dim, ep_dim];
        data_std = permute(sub_data, perm_order);  % [channels x time x episodes]
        
        data_selected.(level) = data_std;
        fprintf('Extracted channels for %s: original size %s, standardized size %s\n', ...
            level, mat2str(dims), mat2str(size(data_selected.(level))));
    end
end
