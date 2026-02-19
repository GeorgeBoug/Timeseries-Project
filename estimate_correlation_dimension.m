function corr_dim = estimate_correlation_dimension(series, embedding_dim, tau)
%% Estimate correlation dimension using Grassberger-Procaccia algorithm
% Inputs:
%   series - time series
%   embedding_dim - embedding dimension
%   tau - time delay (default: 1)
% Output:
%   corr_dim - estimated correlation dimension

    if nargin < 3
        tau = 1;
    end
    
    n = length(series);
    m = n - (embedding_dim - 1) * tau;
    
    if m < 100
        corr_dim = NaN;
        return;
    end
    
    % Create delay embedding
    embedded = zeros(m, embedding_dim);
    for i = 1:embedding_dim
        embedded(:, i) = series((i-1)*tau + 1 : (i-1)*tau + m);
    end
    
    % Compute pairwise distances (sample if too many points)
    max_points = 1000;
    if m > max_points
        indices = randperm(m, max_points);
        embedded = embedded(indices, :);
        m = max_points;
    end
    
    distances = pdist(embedded);
    distances = distances(distances > 0);  % Remove zero distances
    
    % Range of r values (logarithmically spaced)
    r_min = prctile(distances, 1);
    r_max = prctile(distances, 99);
    r_values = logspace(log10(r_min), log10(r_max), 20);
    
    % Compute correlation integral C(r)
    C_r = zeros(size(r_values));
    
    for i = 1:length(r_values)
        r = r_values(i);
        C_r(i) = sum(distances < r) / (length(distances) * (length(distances) - 1) / 2);
    end
    
    % Fit log(C(r)) vs log(r) in linear region
    % Find linear region (avoid saturation and noise)
    log_r = log10(r_values);
    log_C = log10(C_r + eps);  % Add small epsilon to avoid log(0)
    
    % Use middle 60% of points for fitting
    start_idx = round(0.2 * length(log_r));
    end_idx = round(0.8 * length(log_r));
    
    if end_idx > start_idx
        p = polyfit(log_r(start_idx:end_idx), log_C(start_idx:end_idx), 1);
        corr_dim = p(1);  % Slope is the correlation dimension
    else
        corr_dim = NaN;
    end
    
    % Ensure reasonable bounds
    if corr_dim < 0 || corr_dim > embedding_dim * 2
        corr_dim = NaN;
    end
end
