function [values_pre, values_post] = compute_nonlinear_measures_all(...
    segments_pre, segments_post, measure_name, measure_params, sampling_freq)
%% Compute nonlinear measures for all preTMS and postTMS segments
% Inputs:
%   segments_pre - structure with preTMS segments
%   segments_post - structure with postTMS segments
%   measure_name - name of selected nonlinear measure
%   measure_params - parameters for the measure
%   sampling_freq - sampling frequency
% Outputs:
%   values_pre - structure with measure values for preTMS (each intensity level)
%   values_post - structure with measure values for postTMS (each intensity level)

    values_pre = struct();
    values_post = struct();
    intensity_levels = fieldnames(segments_pre);
    
    fprintf('Computing nonlinear measure "%s" for all segments...\n', measure_name);
    
    for i = 1:length(intensity_levels)
        level = intensity_levels{i};
        
        % Process preTMS segments
        pre_data = segments_pre.(level);
        pre_values = compute_nonlinear_measure_for_segments(pre_data, measure_name, ...
            measure_params, sampling_freq);
        values_pre.(level) = pre_values;
        
        % Process postTMS segments
        post_data = segments_post.(level);
        post_values = compute_nonlinear_measure_for_segments(post_data, measure_name, ...
            measure_params, sampling_freq);
        values_post.(level) = post_values;
        
        fprintf('  %s: preTMS = %d values, postTMS = %d values\n', ...
            level, length(pre_values), length(post_values));
    end
end

function values = compute_nonlinear_measure_for_segments(data, measure_name, ...
    measure_params, sampling_freq)
%% Compute nonlinear measure for all segments in a dataset

    % After preprocessing, data is standardized to [channels x time x episodes]
    dims = size(data);
    if length(dims) ~= 3
        error('Expected 3D data');
    end
    n_channels = dims(1);
    n_episodes = dims(3);
    
    values = [];
    
    for ep = 1:n_episodes
        for ch = 1:n_channels
            % Extract time series for this channel/episode
            series = squeeze(data(ch, :, ep));
            
            series = series(:);  % Ensure column vector
            
            % Compute measure based on type
            switch measure_name
                case 'LocalPred_NRMSE_step'
                    % Local nonlinear prediction NRMSE at a fixed horizon step
                    step = measure_params.step;
                    tau = measure_params.tau;
                    m = measure_params.m;
                    nnei = measure_params.nnei;
                    q = measure_params.q;
                    Tmax = measure_params.Tmax;
                    n_fit = measure_params.n_fit;
                    n_forecast = measure_params.n_forecast;
                    
                    if length(series) >= (n_fit + n_forecast)
                        xuse = series(1:(n_fit + n_forecast));
                        nlast = n_forecast;
                    else
                        xuse = series;
                        nlast = min(n_forecast, max(50, round(length(xuse) / 3)));
                    end
                    
                    if step > Tmax
                        Tmax = step;
                    end
                    
                    try
                        nrmseV = localpredictnrmse(xuse, nlast, tau, m, Tmax, nnei, q, []);
                        values = [values; nrmseV(step)];
                    catch
                        values = [values; NaN];
                    end
                    
                otherwise
                    error('Unknown nonlinear measure: %s', measure_name);
            end
        end
    end
end
