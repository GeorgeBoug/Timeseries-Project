function [measures_all, selected_measure, measure_params] = ...
    pilot_nonlinear_analysis(pilot_series, sampling_freq)
%% Pilot nonlinear analysis (uses labnonlinear/)
% Measures on the 6 pilot time series:
% - Mutual Information I(tau) and lag of first minimum
% - False Nearest Neighbors (FNN) and embedding dimension where FNN < 1%
% - Local prediction NRMSE(T) using a local linear model

    %#ok<NASGU>
    n_series = length(pilot_series);

    max_lag_mi = 50;
    max_embedding_dim = 10;
    tau = 1;
    escape = 10;
    theiler = 0;

    % Align nonlinear prediction with the linear train/test split
    n_fit = 600;
    n_forecast = 300; % last part used as test set
    Tmax = 10;
    m = 3;
    nnei = 10;
    q = 3; % q>=m -> local linear OLS

    measures_all = struct();
    measures_all.mutual_info = cell(n_series, 1);           % each: (tmax+1)x2 [lag, MI]
    measures_all.mi_first_min_lag = NaN(n_series, 1);
    measures_all.fnn = cell(n_series, 1);                   % each: mmaxx2 [m, FNN]
    measures_all.fnn_embedding_dim = NaN(n_series, 1);
    measures_all.localpred_nrmse = NaN(n_series, Tmax);     % NRMSE(T), T=1..Tmax

    fprintf('Computing nonlinear measures (labnonlinear)...\n');

    % Ensure folder for figures
    outDir = fullfile(pwd, 'figures', 'pilot_nonlinear');
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    fig_nl = figure('Name', 'Pilot Nonlinear Analysis (labnonlinear)', 'Position', [100, 100, 1400, 900]);
    colors = lines(n_series);

    for i = 1:n_series
        x = pilot_series{i}(:);
        fprintf('Processing series %d (length: %d)...\n', i, length(x));

        % 1) Mutual information curve + first minimum lag
        try
            mutM = mutualinformation(x, max_lag_mi, [], [], 'b'); % no plot (tittxt empty)
            measures_all.mutual_info{i} = mutM;
            miV = mutM(:, 2);

            firstMinLag = NaN;
            for k = 2:(length(miV) - 1) % skip lag 0
                if miV(k) < miV(k-1) && miV(k) < miV(k+1)
                    firstMinLag = mutM(k, 1);
                    break;
                end
            end
            if isnan(firstMinLag)
                [~, idx] = min(miV(2:end));
                firstMinLag = mutM(idx + 1, 1);
            end
            measures_all.mi_first_min_lag(i) = firstMinLag;
        catch ME
            warning('Series %d: mutualinformation failed: %s', i, ME.message);
            measures_all.mutual_info{i} = [];
            measures_all.mi_first_min_lag(i) = NaN;
        end

        % 2) False nearest neighbors
        try
            fnnM = falsenearest(x, tau, max_embedding_dim, escape, theiler, []);
            measures_all.fnn{i} = fnnM;
            fnnV = fnnM(:, 2);
            below = find(fnnV < 0.01, 1);
            if ~isempty(below)
                measures_all.fnn_embedding_dim(i) = fnnM(below, 1);
            else
                measures_all.fnn_embedding_dim(i) = max_embedding_dim;
            end
        catch ME
            warning('Series %d: falsenearest failed: %s', i, ME.message);
            measures_all.fnn{i} = [];
            measures_all.fnn_embedding_dim(i) = NaN;
        end

        % 3) Local prediction NRMSE(T)
        try
            if length(x) >= (n_fit + n_forecast)
                xuse = x(1:(n_fit + n_forecast));
                nlast_use = n_forecast;
            else
                xuse = x;
                nlast_use = min(n_forecast, max(50, round(length(xuse) / 3)));
            end
            nrmseV = localpredictnrmse(xuse, nlast_use, tau, m, Tmax, nnei, q, []);
            measures_all.localpred_nrmse(i, :) = nrmseV(:)';
        catch ME
            warning('Series %d: localpredictnrmse failed: %s', i, ME.message);
            measures_all.localpred_nrmse(i, :) = NaN;
        end
    end

    % --- Plots (comparisons) ---
    subplot(2, 3, 1);
    hold on;
    for i = 1:n_series
        if ~isempty(measures_all.mutual_info{i})
            plot(measures_all.mutual_info{i}(:, 1), measures_all.mutual_info{i}(:, 2), ...
                'Color', colors(i, :), 'LineWidth', 1.5, 'DisplayName', sprintf('Series %d', i));
        end
    end
    xlabel('Lag');
    ylabel('I(\tau)');
    title('Mutual Information');
    legend('show');
    grid on;
    hold off;

    subplot(2, 3, 2);
    hold on;
    for i = 1:n_series
        if ~isempty(measures_all.fnn{i})
            plot(measures_all.fnn{i}(:, 1), measures_all.fnn{i}(:, 2), ...
                'Color', colors(i, :), 'LineWidth', 1.5, 'Marker', 'o', 'DisplayName', sprintf('Series %d', i));
        end
    end
    plot([1, max_embedding_dim], [0.01, 0.01], 'r--', 'LineWidth', 1);
    xlabel('Embedding dimension m');
    ylabel('FNN(m)');
    title('False Nearest Neighbors');
    legend('show');
    grid on;
    hold off;

    subplot(2, 3, 3);
    hold on;
    for i = 1:n_series
        plot(1:Tmax, measures_all.localpred_nrmse(i, :), ...
            'Color', colors(i, :), 'LineWidth', 1.5, 'Marker', 'o', 'DisplayName', sprintf('Series %d', i));
    end
    xlabel('Prediction horizon T');
    ylabel('NRMSE(T)');
    title('Local prediction NRMSE(T)');
    legend('show');
    grid on;
    hold off;

    subplot(2, 3, 4);
    bar(measures_all.mi_first_min_lag);
    xlabel('Series');
    ylabel('Lag');
    title('MI first minimum lag');
    grid on;

    subplot(2, 3, 5);
    bar(measures_all.fnn_embedding_dim);
    xlabel('Series');
    ylabel('m');
    title('Embedding dim where FNN < 1%');
    grid on;

    subplot(2, 3, 6);
    step = 5;
    bar(measures_all.localpred_nrmse(:, step));
    xlabel('Series');
    ylabel(sprintf('NRMSE(T=%d)', step));
    title('Local prediction NRMSE at step 5');
    grid on;

    % Save pilot nonlinear overview
    try
        saveas(fig_nl, fullfile(outDir, 'pilot_nonlinear_overview.png'));
    catch
    end

    % Default selection for full analysis
    selected_measure = 'LocalPred_NRMSE_step';
    measure_params.step = 5;
    measure_params.tau = tau;
    measure_params.m = m;
    measure_params.nnei = nnei;
    measure_params.q = q;
    measure_params.Tmax = Tmax;
    measure_params.n_fit = n_fit;
    measure_params.n_forecast = n_forecast;

    fprintf('Selected nonlinear measure: %s (step=%d)\n', selected_measure, measure_params.step);
end
