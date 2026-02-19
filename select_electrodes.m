function [electrode_indices, electrode_names] = select_electrodes(aem_numbers, channel_names)
%% Select two electrodes based on AEM numbers using mod(x,60)+1 rule
% Inputs:
%   aem_numbers - array of AEM numbers (e.g., [1, 1])
%   channel_names - cell array of channel names
% Outputs:
%   electrode_indices - array of two electrode indices (1-based)
%   electrode_names - cell array of two electrode names

    n_channels = length(channel_names);
    
    % Calculate electrode indices using mod(x,60)+1
    idx1 = mod(aem_numbers(1), 60) + 1;
    idx2 = mod(aem_numbers(2), 60) + 1;
    
    % Adjust if indices exceed number of channels
    if idx1 > n_channels
        idx1 = mod(idx1 - 1, n_channels) + 1;
    end
    if idx2 > n_channels
        idx2 = mod(idx2 - 1, n_channels) + 1;
    end
    
    % Adjust if two indices coincide
    if idx1 == idx2
        idx2 = mod(idx2, n_channels) + 1;
        if idx2 == idx1
            idx2 = mod(idx2, n_channels) + 1;
        end
    end
    
    electrode_indices = [idx1, idx2];
    electrode_names = {channel_names{idx1}, channel_names{idx2}};
    
    fprintf('Electrode selection:\n');
    fprintf('  AEM number 1: %d -> electrode index %d (%s)\n', aem_numbers(1), idx1, electrode_names{1});
    fprintf('  AEM number 2: %d -> electrode index %d (%s)\n', aem_numbers(2), idx2, electrode_names{2});
end
