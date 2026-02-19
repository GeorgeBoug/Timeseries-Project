function rejected_episodes = visualize_and_reject_episodes(data_selected, intensity_levels, sampling_freq)
%% Visualize all episodes and allow manual rejection of artifacts/epileptiform discharges
% Inputs:
%   data_selected - structure with selected channels for each intensity level
%   intensity_levels - cell array of intensity level names
%   sampling_freq - sampling frequency in Hz
% Outputs:
%   rejected_episodes - structure with rejected episode indices for each intensity level

    rejected_episodes = struct();
    
    for i = 1:length(intensity_levels)
        level = intensity_levels{i};
        data = data_selected.(level);
        
        % After extraction, data is standardized to [channels x time x episodes]
        dims = size(data);
        if length(dims) ~= 3
            error('Expected 3D data for %s', level);
        end
        n_channels = dims(1);
        n_samples  = dims(2);
        n_episodes = dims(3);
        
        % Reshape to [episodes x channels x time] for easier iteration
        data_reshaped = permute(data, [3, 1, 2]);
        
        fprintf('\nVisualizing episodes for intensity level: %s\n', level);
        fprintf('Total episodes: %d, Channels: %d, Samples: %d\n', n_episodes, n_channels, n_samples);
        
        % Create figure for visualization
        fig = figure('Name', sprintf('Episode Review: %s', level), ...
            'Position', [100, 100, 1200, 800]);
        
        rejected = false(n_episodes, 1);
        current_episode = 1;
        
        while current_episode <= n_episodes
            % Plot current episode for both channels
            clf(fig);
            
            for ch = 1:n_channels
                subplot(n_channels, 1, ch);
                episode_data = squeeze(data_reshaped(current_episode, ch, :));
                time_axis = (0:length(episode_data)-1) / sampling_freq * 1000;  % Convert to ms
                
                plot(time_axis, episode_data, 'b-', 'LineWidth', 1);
                title(sprintf('Episode %d/%d - Channel %d - %s', ...
                    current_episode, n_episodes, ch, level));
                xlabel('Time (ms)');
                ylabel('Amplitude');
                grid on;
                
                if rejected(current_episode)
                    title(sprintf('Episode %d/%d - Channel %d - %s [REJECTED]', ...
                        current_episode, n_episodes, ch, level), 'Color', 'r');
                end
            end
            
            % Remove any previous instruction box
            delete(findall(fig, 'Type', 'textboxshape'));

            % Instructions (placed on the right side so they don't cover signals)
            msg = sprintf(['Episode %d/%d\n\n' ...
                'Controls:\n' ...
                '  [Space] or\n' ...
                ' [Right Arrow]-\n ' ...
                ' Next episode\n' ...
                '  \n [Left Arrow] -\n' ...
                ' Previous episode\n' ...
                '  \n [R] - Toggle reject/\n' ...
                'accept\n' ...
                '  \n [Q] - Quit and save\n' ...
                '  \n [A] - Accept all\n' ...
                'remaining\n' ...
                '\n \nRejected: %d/%d'], ...
                current_episode, n_episodes, sum(rejected), n_episodes);
            
            % [x y w h] in normalized units; put on the right side
            annotation('textbox', [0.91, 0.20, 0.14, 0.33], ...
                'String', msg, 'FontSize', 7, ...
                'BackgroundColor', 'w', 'EdgeColor', 'k');
            
            % Wait for user key press
            waitforbuttonpress;
            key = get(fig, 'CurrentKey');  % works for letters and arrow keys
            if isempty(key)
                continue;
            end
            
            switch key
                case {'space', 'rightarrow'}
                    current_episode = current_episode + 1;
                case 'leftarrow'
                    current_episode = max(1, current_episode - 1);
                case {'r', 'R'}
                    rejected(current_episode) = ~rejected(current_episode);
                case {'q', 'Q'}
                    break;
                case {'a', 'A'}
                    % Accept all remaining: just exit loop without adding more rejections
                    break;
            end
        end
        
        rejected_episodes.(level) = find(rejected);
        fprintf('Rejected %d episodes out of %d for %s\n', ...
            sum(rejected), n_episodes, level);
        
        close(fig);
    end
    
    fprintf('\nManual rejection complete.\n');
end
