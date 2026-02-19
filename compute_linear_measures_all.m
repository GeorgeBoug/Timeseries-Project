function [values_pre, values_post] = compute_linear_measures_all(...
    segments_pre, segments_post, measure_name, measure_params, ...
    sampling_freq, n_fit, n_forecast, forecast_horizons)
%% Compute linear measures for all preTMS and postTMS segments
% Inputs:
%   segments_pre - structure with preTMS segments
%   segments_post - structure with postTMS segments
%   measure_name - name of selected linear measure
%   measure_params - parameters for the measure
%   sampling_freq - sampling frequency
%   n_fit, n_forecast, forecast_horizons - forecasting parameters
% Outputs:
%   values_pre - structure with measure values for preTMS (each intensity level)
%   values_post - structure with measure values for postTMS (each intensity level)

    values_pre = struct();
    values_post = struct();
    intensity_levels = fieldnames(segments_pre);
    
    fprintf('Computing linear measure "%s" for all segments...\n', measure_name);
    
    for i = 1:length(intensity_levels)
        level = intensity_levels{i};
        
        % Process preTMS segments
        pre_data = segments_pre.(level);
        pre_values = compute_measure_for_segments(pre_data, measure_name, ...
            measure_params, sampling_freq, n_fit, n_forecast, forecast_horizons);
        values_pre.(level) = pre_values;
        
        % Process postTMS segments
        post_data = segments_post.(level);
        post_values = compute_measure_for_segments(post_data, measure_name, ...
            measure_params, sampling_freq, n_fit, n_forecast, forecast_horizons);
        values_post.(level) = post_values;
        
        fprintf('  %s: preTMS = %d values, postTMS = %d values\n', ...
            level, length(pre_values), length(post_values));
    end
end

function values = compute_measure_for_segments(data, measure_name, ...
    measure_params, sampling_freq, n_fit, n_forecast, forecast_horizons)
%% Compute measure for all segments in a dataset

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
            series = series(:);
            
            % Compute measure based on type
            switch measure_name
                case 'NRMSE_step5'
                    % ARMA prediction NRMSE at a fixed horizon using lablinear/predictARMAnrmse
                    horizon = measure_params.forecast_horizon;
                    if length(series) >= (n_fit + n_forecast)
                        try
                            xuse = series(1:(n_fit+n_forecast));
                            Tmax = max(horizon, 10);
                            nlast = n_forecast;
                            
                            % Heuristic p,q from ACF/PACF (lablinear)
                            max_lag = min(20, floor(length(xuse) / 4));
                            acM = autocorrelation(xuse, max_lag, [], 'b');
                            pautV = parautocor(xuse, max_lag);
                            pacf_threshold = 0.2;
                            acf_threshold  = 0.2;
                            
                            p = 0;
                            for lag = 1:min(10, length(pautV))
                                if abs(pautV(lag)) < pacf_threshold
                                    p = max(lag-1, 0);
                                    break;
                                end
                            end
                            q = 0;
                            for lag = 2:min(11, size(acM,1)) % acM includes lag 0
                                if abs(acM(lag,2)) < acf_threshold
                                    q = max((lag-1)-1, 0);
                                    break;
                                end
                            end
                            
                            nrmseV = predictARMAnrmse(xuse, p, q, Tmax, nlast, []);
                            values = [values; nrmseV(horizon)];
                        catch
                            values = [values; NaN];
                        end
                    else
                        values = [values; NaN];
                    end
                    
                case 'ACF_lag10'
                    lag = measure_params.lag;
                    if length(series) > lag
                        try
                            acM = autocorrelation(series, lag, [], 'b');
                            values = [values; acM(lag+1,2)]; % includes lag 0 at row 1
                        catch
                            values = [values; NaN];
                        end
                    else
                        values = [values; NaN];
                    end
                    
                case 'PACF_lag5'
                    lag = measure_params.lag;
                    if length(series) > lag
                        try
                            pautV = parautocor(series, lag);
                            values = [values; pautV(lag)];
                        catch
                            values = [values; NaN];
                        end
                    else
                        values = [values; NaN];
                    end
                    
                otherwise
                    error('Unknown linear measure: %s', measure_name);
            end
        end
    end
end
