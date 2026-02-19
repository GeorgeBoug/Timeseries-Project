function nrmse = compute_nonlinear_forecasting_nrmse(series, embedding_dim, tau, n_fit, n_forecast)
%% Compute NRMSE for local linear nonlinear forecasting
% Uses delay embedding and local linear model for forecasting

    n = length(series);
    if n_fit + n_forecast > n
        n_forecast = n - n_fit;
    end
    
    % Create delay embedding for training data
    m_fit = n_fit - (embedding_dim - 1) * tau;
    embedded_fit = zeros(m_fit, embedding_dim);
    for i = 1:embedding_dim
        embedded_fit(:, i) = series((i-1)*tau + 1 : (i-1)*tau + m_fit);
    end
    targets_fit = series(embedding_dim * tau + 1 : embedding_dim * tau + m_fit);
    
    % Forecast
    forecasts = zeros(n_forecast, 1);
    k_neighbors = min(embedding_dim + 1, floor(m_fit / 2));  % Number of neighbors
    
    for t = 1:n_forecast
        % Current state vector
        current_state = zeros(1, embedding_dim);
        for i = 1:embedding_dim
            idx = n_fit + t - (embedding_dim - i) * tau;
            if idx > 0 && idx <= n
                current_state(i) = series(idx);
            else
                current_state(i) = 0;
            end
        end
        
        % Find k nearest neighbors in training set
        distances = pdist2(embedded_fit, current_state);
        [~, nn_indices] = mink(distances, k_neighbors);
        
        % Local linear model: fit linear regression to neighbors
        X_neighbors = [ones(k_neighbors, 1), embedded_fit(nn_indices, :)];
        y_neighbors = targets_fit(nn_indices);
        
        % Weighted least squares (weights based on distance)
        weights = 1 ./ (distances(nn_indices) + eps);
        weights = weights / sum(weights);
        W = diag(weights);
        
        % Solve weighted least squares
        beta = (X_neighbors' * W * X_neighbors) \ (X_neighbors' * W * y_neighbors);
        
        % Forecast
        forecasts(t) = beta(1) + current_state * beta(2:end);
    end
    
    % True values
    true_values = series(n_fit + 1 : n_fit + n_forecast);
    
    % Compute NRMSE
    rmse = sqrt(mean((true_values - forecasts).^2));
    nrmse = rmse / (max(true_values) - min(true_values));
end
