function [acf_all, pacf_all, nrmse_all, selected_measure, measure_params] = ...
    pilot_linear_analysis(pilot_series, sampling_freq, n_fit, n_forecast, forecast_horizons)
%% Pilot linear analysis (LABLINEAR): ACF, PACF, ARMA modeling, forecasting
% Inputs:
%   pilot_series - cell array of 6 time series
%   sampling_freq - sampling frequency in Hz
%   n_fit - number of samples for model fitting
%   n_forecast - number of samples for forecasting
%   forecast_horizons - array of forecast step-ahead values
% Outputs:
%   acf_all - cell array of ACF values for each series
%   pacf_all - cell array of PACF values for each series
%   nrmse_all - matrix of NRMSE values [series x horizon]
%   selected_measure - name of selected linear measure
%   measure_params - parameters for selected measure

    %#ok<NASGU>
    n_series = length(pilot_series);
    max_lag = min(100, floor(length(pilot_series{1}) / 4));
    Tmax = max(forecast_horizons);
    nlast = n_forecast; % test set size for predictARMAnrmse
    
    acf_all = cell(n_series, 1);   % each: (tmax+1)x2 [lag, acf]
    pacf_all = cell(n_series, 1);  % each: (tmax)x2 [lag(1..tmax), pacf]
    nrmse_all = NaN(n_series, length(forecast_horizons));
    
    fprintf('Computing ACF/PACF and fitting ARMA models (lablinear)...\n');
    
    % Ensure folder for figures
    outDir = fullfile(pwd, 'figures', 'pilot_linear');
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    % ACF/PACF per series
    fig_acf = figure('Name', 'Pilot ACF/PACF Analysis (lablinear)', 'Position', [100, 100, 1400, 900]);
    for i = 1:n_series
        x = pilot_series{i}(:);
        
        % ACF (no plot from lab function when tittxt empty)
        acM = autocorrelation(x, max_lag, [], 'b');
        acf_all{i} = acM;
        
        % PACF
        pautV = parautocor(x, max_lag);
        pacf_all{i} = [(1:max_lag)', pautV(:)];
        
        subplot(3, 2, i);
        yyaxis left
        stem(acM(:,1), acM(:,2), 'b', 'MarkerSize', 4);
        ylabel('ACF');
        yyaxis right
        stem(1:max_lag, pautV, 'g', 'MarkerSize', 4);
        ylabel('PACF');
        xlabel('Lag');
        title(sprintf('Series %d: ACF (blue) / PACF (green)', i));
        grid on;
    end
    % Save the per-series ACF/PACF overview
    try
        saveas(fig_acf, fullfile(outDir, 'pilot_acf_pacf_per_series.png'));
    catch
    end
    
    % Compare ACF, PACF, and NRMSE curves
    fprintf('\nComparing measures across series...\n');
    
    % Plot comparison of ACF curves
    fig_acf_comp = figure('Name', 'ACF Comparison', 'Position', [200, 200, 1000, 600]);
    hold on;
    colors = lines(n_series);
    for i = 1:n_series
        plot(acf_all{i}(:,1), acf_all{i}(:,2), 'Color', colors(i,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Series %d', i));
    end
    xlim([0 max_lag]);
    plot([0, max_lag], [0.2, 0.2], 'r--', 'LineWidth', 1);
    plot([0, max_lag], [-0.2, -0.2], 'r--', 'LineWidth', 1);
    xlabel('Lag');
    ylabel('ACF');
    title('ACF Comparison Across Pilot Series');
    legend('show');
    grid on;
    hold off;
    % Save ACF comparison
    try
        saveas(fig_acf_comp, fullfile(outDir, 'acf_comparison.png'));
    catch
    end
    
    % Plot comparison of PACF curves
    fig_pacf_comp = figure('Name', 'PACF Comparison', 'Position', [300, 300, 1000, 600]);
    hold on;
    for i = 1:n_series
        plot(pacf_all{i}(:,1), pacf_all{i}(:,2), 'Color', colors(i,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Series %d', i));
    end
    xlim([1 max_lag]);
    plot([1, max_lag], [0.2, 0.2], 'r--', 'LineWidth', 1);
    plot([1, max_lag], [-0.2, -0.2], 'r--', 'LineWidth', 1);
    xlabel('Lag');
    ylabel('PACF');
    title('PACF Comparison Across Pilot Series');
    legend('show');
    grid on;
    hold off;
    % Save PACF comparison
    try
        saveas(fig_pacf_comp, fullfile(outDir, 'pacf_comparison.png'));
    catch
    end
    
    % --- Fit ARMA and compute prediction NRMSE curves using lablinear ---
    fprintf('\nARMA prediction error (lablinear/predictARMAnrmse) on training=%d, test=%d...\n', n_fit, n_forecast);
    for i = 1:n_series
        x = pilot_series{i}(:);
        if length(x) < (n_fit + n_forecast)
            warning('Series %d is too short (%d) for %d+%d split. Skipping NRMSE.', i, length(x), n_fit, n_forecast);
            continue;
        end
        xuse = x(1:(n_fit+n_forecast));
        
        % Heuristic p,q from early ACF/PACF
        pacf_threshold = 0.2;
        acf_threshold  = 0.2;
        pautV = pacf_all{i}(:,2);
        acfV  = acf_all{i}(:,2);
        
        p = 0;
        for lag = 1:min(10, length(pautV))
            if abs(pautV(lag)) < pacf_threshold
                p = max(lag-1, 0);
                break;
            end
        end
        q = 0;
        for lag = 2:min(11, length(acfV)) % acfV includes lag 0 at index 1
            if abs(acfV(lag)) < acf_threshold
                q = max((lag-1)-1, 0); % map index to lag then -1
                break;
            end
        end
        
        % Residual whiteness check on training set (Ljung-Box via lablinear)
        try
            [~,~,~,~,~,~,~,xpreM] = fitARMA(xuse(1:n_fit), p, q, 1);
            res = xuse(1:n_fit) - xpreM(:,1);
            [hV, pV] = portmanteauLB(res, 10, 0.05, []); %#ok<ASGLU>
            fprintf('Series %d: ARMA(%d,%d), Portmanteau p(lag=10)=%.4f (white: %d)\n', ...
                i, p, q, pV(end), hV(end)==0);
        catch ME
            warning('Series %d: residual test failed: %s', i, ME.message);
        end
        
        % Prediction NRMSE(T), T=1..Tmax on last nlast samples
        try
            [nrmseV] = predictARMAnrmse(xuse, p, q, Tmax, nlast, []); % no plot
            for hi = 1:length(forecast_horizons)
                T = forecast_horizons(hi);
                nrmse_all(i, hi) = nrmseV(T);
            end
        catch ME
            warning('Series %d: predictARMAnrmse failed: %s', i, ME.message);
        end
    end
    
    % Plot comparison of NRMSE curves
    fig_nrmse_comp = figure('Name', 'NRMSE Comparison', 'Position', [400, 400, 1000, 600]);
    hold on;
    for i = 1:n_series
        plot(forecast_horizons, nrmse_all(i, :), 'Color', colors(i,:), ...
            'LineWidth', 1.5, 'Marker', 'o', 'DisplayName', sprintf('Series %d', i));
    end
    xlabel('Forecast Horizon (steps)');
    ylabel('NRMSE');
    title('NRMSE Comparison Across Pilot Series');
    legend('show');
    grid on;
    hold off;
    % Save NRMSE comparison
    try
        saveas(fig_nrmse_comp, fullfile(outDir, 'nrmse_comparison.png'));
    catch
    end
    
    % Select one linear measure (user can modify this selection logic)
    % Default: Use NRMSE at step 5
    selected_measure = 'NRMSE_step5';
    measure_params.forecast_horizon = 5;
    measure_params.n_fit = n_fit;
    measure_params.n_forecast = n_forecast;
    measure_params.p = []; % let full-analysis pick p,q per-series by same heuristic
    measure_params.q = [];
    
    fprintf('Selected linear measure: %s (forecast horizon = %d)\n', ...
        selected_measure, measure_params.forecast_horizon);
    
    % Alternative selections (commented out):
    % selected_measure = 'ACF_lag10';
    % measure_params.lag = 10;
    %
    % selected_measure = 'PACF_lag5';
    % measure_params.lag = 5;
end
