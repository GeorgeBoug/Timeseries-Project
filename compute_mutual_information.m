function mi = compute_mutual_information(series, lag)
%% Compute mutual information between X(t) and X(t+lag)
% Uses histogram-based estimation with adaptive binning

    n = length(series);
    if lag >= n
        mi = 0;
        return;
    end
    
    x = series(1:n-lag);
    y = series(1+lag:n);
    
    % Adaptive binning using Freedman-Diaconis rule
    n_bins_x = ceil((max(x) - min(x)) / (2 * iqr(x) / (n-lag)^(1/3)));
    n_bins_y = ceil((max(y) - min(y)) / (2 * iqr(y) / (n-lag)^(1/3)));
    
    n_bins_x = max(5, min(20, n_bins_x));
    n_bins_y = max(5, min(20, n_bins_y));
    
    % Compute joint and marginal histograms
    [counts_xy, ~, ~] = histcounts2(x, y, n_bins_x, n_bins_y);
    [counts_x, ~] = histcounts(x, n_bins_x);
    [counts_y, ~] = histcounts(y, n_bins_y);
    
    % Normalize to probabilities
    p_xy = counts_xy / sum(counts_xy(:));
    p_x = counts_x / sum(counts_x);
    p_y = counts_y / sum(counts_y);
    
    % Compute mutual information: MI = sum(p(x,y) * log(p(x,y) / (p(x)*p(y))))
    mi = 0;
    for i = 1:size(p_xy, 1)
        for j = 1:size(p_xy, 2)
            if p_xy(i,j) > 0 && p_x(i) > 0 && p_y(j) > 0
                mi = mi + p_xy(i,j) * log2(p_xy(i,j) / (p_x(i) * p_y(j)));
            end
        end
    end
end
