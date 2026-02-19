function [data_all, channel_names] = load_eeg_data(mat_file_paths, channels_file, intensity_levels)
%% Load EEG data from .mat files and channel names from text file
% Inputs:
%   mat_file_paths - cell array of full paths to .mat files (one per intensity level)
%   channels_file - path to Channels.txt file
%   intensity_levels - cell array of intensity level names
% Outputs:
%   data_all - structure with fields for each intensity level containing cMF arrays
%   channel_names - cell array of channel names (excluding first 4 non-EEG channels)

    % Load channel names
    fprintf('Loading channel names from %s...\n', channels_file);
    fid = fopen(channels_file, 'r');
    if fid == -1
        error('Cannot open file: %s', channels_file);
    end
    
    channel_names_raw = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    channel_names_raw = channel_names_raw{1};
    
    % Remove quotes if present and trim whitespace
    for i = 1:length(channel_names_raw)
        channel_names_raw{i} = strtrim(channel_names_raw{i});
        % Remove surrounding single or double quotes if present
        if length(channel_names_raw{i}) >= 2
            if (channel_names_raw{i}(1) == '''' && channel_names_raw{i}(end) == '''') || ...
               (channel_names_raw{i}(1) == '"' && channel_names_raw{i}(end) == '"')
                channel_names_raw{i} = channel_names_raw{i}(2:end-1);
            end
        end
    end
    
    % Skip first 4 channels (non-EEG)
    channel_names = channel_names_raw(5:end);
    fprintf('Loaded %d EEG channel names (skipped first 4 non-EEG channels)\n', length(channel_names));
    
    % Load .mat files for each intensity level
    data_all = struct();
    
    for i = 1:length(intensity_levels)
        level = intensity_levels{i};
        
        if i <= length(mat_file_paths)
            mat_file = mat_file_paths{i};
        else
            error('Not enough .mat file paths provided for intensity level: %s', level);
        end
        
        if ~exist(mat_file, 'file')
            error('Cannot find .mat file: %s', mat_file);
        end
        
        fprintf('Loading %s...\n', mat_file);
        loaded_data = load(mat_file);
        
        % Extract cMF array
        if isfield(loaded_data, 'cMF')
            data_all.(level) = loaded_data.cMF;
        elseif isfield(loaded_data, 'cMF_array')
            data_all.(level) = loaded_data.cMF_array;
        else
            % Try to find any variable containing 'cMF' in its name
            vars = fieldnames(loaded_data);
            cMF_idx = contains(vars, 'cMF', 'IgnoreCase', true);
            if any(cMF_idx)
                data_all.(level) = loaded_data.(vars{find(cMF_idx, 1)});
            else
                % List available variables to help debug
                fprintf('Available variables in file: %s\n', strjoin(vars, ', '));
                error('Cannot find cMF array in file: %s', mat_file);
            end
        end
        
        fprintf('Loaded %s: size %s\n', level, mat2str(size(data_all.(level))));
    end
end
