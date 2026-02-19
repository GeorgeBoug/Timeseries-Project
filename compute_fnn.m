function fnn_ratio = compute_fnn(series, embedding_dim, tau)
%% Compute False Nearest Neighbors ratio for given embedding dimension
% Uses the method by Kennel et al. (1992)

    n = length(series);
    if embedding_dim * tau >= n
        fnn_ratio = 1;
        return;
    end
    
    % Create delay embedding
    m = n - (embedding_dim - 1) * tau;
    embedded = zeros(m, embedding_dim);
    for i = 1:embedding_dim
        embedded(:, i) = series((i-1)*tau + 1 : (i-1)*tau + m);
    end
    
    if embedding_dim == 1
        % For dimension 1, compare with dimension 2
        embedded_next = zeros(m - tau, embedding_dim + 1);
        for i = 1:embedding_dim + 1
            embedded_next(:, i) = series((i-1)*tau + 1 : (i-1)*tau + m - tau);
        end
        
        % Find nearest neighbors in dimension 1
        distances = pdist2(embedded(1:end-tau, :), embedded(1:end-tau, :));
        distances(logical(eye(size(distances)))) = inf;  % Exclude self
        
        [~, nn_indices] = min(distances, [], 2);
        
        % Check if they are still neighbors in dimension 2
        fnn_count = 0;
        for i = 1:length(nn_indices)
            dist_dim1 = norm(embedded(i, :) - embedded(nn_indices(i), :));
            dist_dim2 = norm(embedded_next(i, :) - embedded_next(nn_indices(i), :));
            
            if dist_dim1 > 0
                ratio = abs(dist_dim2 - dist_dim1) / dist_dim1;
                if ratio > 15 || dist_dim2 / dist_dim1 > 2
                    fnn_count = fnn_count + 1;
                end
            end
        end
        
        fnn_ratio = fnn_count / length(nn_indices);
    else
        % For higher dimensions, compare with next dimension
        embedded_next = zeros(m - tau, embedding_dim + 1);
        for i = 1:embedding_dim + 1
            embedded_next(:, i) = series((i-1)*tau + 1 : (i-1)*tau + m - tau);
        end
        
        % Find nearest neighbors in current dimension
        distances = pdist2(embedded(1:end-tau, :), embedded(1:end-tau, :));
        distances(logical(eye(size(distances)))) = inf;
        
        [~, nn_indices] = min(distances, [], 2);
        
        % Check if they are still neighbors in next dimension
        fnn_count = 0;
        for i = 1:length(nn_indices)
            dist_dim_d = norm(embedded(i, :) - embedded(nn_indices(i), :));
            dist_dim_d1 = norm(embedded_next(i, :) - embedded_next(nn_indices(i), :));
            
            if dist_dim_d > 0
                ratio = abs(dist_dim_d1 - dist_dim_d) / dist_dim_d;
                if ratio > 15 || dist_dim_d1 / dist_dim_d > 2
                    fnn_count = fnn_count + 1;
                end
            end
        end
        
        fnn_ratio = fnn_count / length(nn_indices);
    end
end
