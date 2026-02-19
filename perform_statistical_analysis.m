function perform_statistical_analysis(linear_values_pre, linear_values_post, ...
    nonlinear_values_pre, nonlinear_values_post, intensity_levels, ...
    linear_measure_name, nonlinear_measure_name)
%% Perform statistical comparisons and create visualizations
% Inputs:
%   linear_values_pre/post - structures with linear measure values
%   nonlinear_values_pre/post - structures with nonlinear measure values
%   intensity_levels - cell array of intensity level names
%   linear_measure_name - name of linear measure
%   nonlinear_measure_name - name of nonlinear measure

    fprintf('\nPerforming statistical analysis...\n');
    
    % Ensure figure output folders
    outDirBase = fullfile(pwd, 'figures', 'statistics');
    if ~exist(outDirBase, 'dir')
        mkdir(outDirBase);
    end
    
    % Process linear measures
    fprintf('\n=== Linear Measure: %s ===\n', linear_measure_name);
    compare_measures(linear_values_pre, linear_values_post, intensity_levels, ...
        linear_measure_name, 'Linear', fullfile(outDirBase, 'linear'));
    
    % Process nonlinear measures
    fprintf('\n=== Nonlinear Measure: %s ===\n', nonlinear_measure_name);
    compare_measures(nonlinear_values_pre, nonlinear_values_post, intensity_levels, ...
        nonlinear_measure_name, 'Nonlinear', fullfile(outDirBase, 'nonlinear'));
end

function compare_measures(values_pre, values_post, intensity_levels, measure_name, measure_type, outDir)
%% Compare preTMS vs postTMS and across intensity levels

    n_levels = length(intensity_levels);
    
    % Extract values into matrices
    pre_matrix = [];
    post_matrix = [];
    level_labels = {};
    
    for i = 1:n_levels
        level = intensity_levels{i};
        pre_vals = values_pre.(level);
        post_vals = values_post.(level);
        
        % Remove NaN values
        pre_vals = pre_vals(~isnan(pre_vals));
        post_vals = post_vals(~isnan(post_vals));
        
        pre_matrix = [pre_matrix; pre_vals];
        post_matrix = [post_matrix; post_vals];
        
        level_labels = [level_labels; repmat({level}, length(pre_vals), 1)];
    end
    
    % Ensure subfolder for this measure type
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    
    % Create figure for all comparisons
    fig = figure('Name', sprintf('%s Measure: %s', measure_type, measure_name), ...
        'Position', [100, 100, 1600, 1200]);
    
    % 1. Compare preTMS vs postTMS within each intensity level
    subplot(3, 3, 1);
    comparison_data = [];
    group_labels = {};
    for i = 1:n_levels
        level = intensity_levels{i};
        pre_vals = values_pre.(level);
        post_vals = values_post.(level);
        pre_vals = pre_vals(~isnan(pre_vals));
        post_vals = post_vals(~isnan(post_vals));
        
        comparison_data = [comparison_data; pre_vals; post_vals];
        group_labels = [group_labels; ...
            repmat({sprintf('%s_pre', level)}, length(pre_vals), 1); ...
            repmat({sprintf('%s_post', level)}, length(post_vals), 1)];
    end
    boxplot(comparison_data, group_labels);
    ylabel(measure_name);
    title('preTMS vs postTMS within each intensity level');
    xtickangle(45);
    grid on;
    
    % Histogram
    subplot(3, 3, 2);
    hold on;
    colors = lines(n_levels);
    for i = 1:n_levels
        level = intensity_levels{i};
        pre_vals = values_pre.(level);
        post_vals = values_post.(level);
        pre_vals = pre_vals(~isnan(pre_vals));
        post_vals = post_vals(~isnan(post_vals));
        
        histogram(pre_vals, 'FaceColor', colors(i,:), 'FaceAlpha', 0.5, ...
            'DisplayName', sprintf('%s_pre', level));
        histogram(post_vals, 'FaceColor', colors(i,:), 'FaceAlpha', 0.3, ...
            'LineStyle', '--', 'DisplayName', sprintf('%s_post', level));
    end
    xlabel(measure_name);
    ylabel('Frequency');
    title('Histogram: preTMS vs postTMS');
    legend('show');
    grid on;
    hold off;
    
    % Scatter plot
    subplot(3, 3, 3);
    hold on;
    for i = 1:n_levels
        level = intensity_levels{i};
        pre_vals = values_pre.(level);
        post_vals = values_post.(level);
        pre_vals = pre_vals(~isnan(pre_vals));
        post_vals = post_vals(~isnan(post_vals));
        
        % Match lengths if needed
        min_len = min(length(pre_vals), length(post_vals));
        scatter(pre_vals(1:min_len), post_vals(1:min_len), ...
            'filled', 'DisplayName', level);
    end
    xlabel('preTMS');
    ylabel('postTMS');
    title('Scatter: preTMS vs postTMS');
    legend('show');
    plot([min([pre_matrix; post_matrix]), max([pre_matrix; post_matrix])], ...
        [min([pre_matrix; post_matrix]), max([pre_matrix; post_matrix])], ...
        'k--', 'LineWidth', 1);
    grid on;
    hold off;
    
    % 2. Compare three intensity levels for postTMS alone
    subplot(3, 3, 4);
    post_data = [];
    post_labels = {};
    for i = 1:n_levels
        level = intensity_levels{i};
        post_vals = values_post.(level);
        post_vals = post_vals(~isnan(post_vals));
        post_data = [post_data; post_vals];
        post_labels = [post_labels; repmat({level}, length(post_vals), 1)];
    end
    boxplot(post_data, post_labels);
    ylabel(measure_name);
    title('postTMS: Comparison across intensity levels');
    grid on;
    
    subplot(3, 3, 5);
    hold on;
    for i = 1:n_levels
        level = intensity_levels{i};
        post_vals = values_post.(level);
        post_vals = post_vals(~isnan(post_vals));
        histogram(post_vals, 'FaceColor', colors(i,:), 'FaceAlpha', 0.6, ...
            'DisplayName', level);
    end
    xlabel(measure_name);
    ylabel('Frequency');
    title('Histogram: postTMS across intensity levels');
    legend('show');
    grid on;
    hold off;
    
    subplot(3, 3, 6);
    scatter_plot_multiple_groups(values_post, intensity_levels, 'postTMS', measure_name);
    
    % 3. Compare three intensity levels for preTMS alone
    subplot(3, 3, 7);
    pre_data = [];
    pre_labels = {};
    for i = 1:n_levels
        level = intensity_levels{i};
        pre_vals = values_pre.(level);
        pre_vals = pre_vals(~isnan(pre_vals));
        pre_data = [pre_data; pre_vals];
        pre_labels = [pre_labels; repmat({level}, length(pre_vals), 1)];
    end
    boxplot(pre_data, pre_labels);
    ylabel(measure_name);
    title('preTMS: Comparison across intensity levels');
    grid on;
    
    subplot(3, 3, 8);
    hold on;
    for i = 1:n_levels
        level = intensity_levels{i};
        pre_vals = values_pre.(level);
        pre_vals = pre_vals(~isnan(pre_vals));
        histogram(pre_vals, 'FaceColor', colors(i,:), 'FaceAlpha', 0.6, ...
            'DisplayName', level);
    end
    xlabel(measure_name);
    ylabel('Frequency');
    title('Histogram: preTMS across intensity levels');
    legend('show');
    grid on;
    hold off;
    
    subplot(3, 3, 9);
    scatter_plot_multiple_groups(values_pre, intensity_levels, 'preTMS', measure_name);
    
    % Save first figure (levels & pre/post)
    try
        fname1 = sprintf('%s_%s_levels_pre_post.png', measure_type, measure_name);
        saveas(fig, fullfile(outDir, fname1));
    catch
    end
    
    % 4. Compare intensity levels using difference (postTMS - preTMS)
    fig2 = figure('Name', sprintf('%s Measure Difference: %s', measure_type, measure_name), ...
        'Position', [200, 200, 1400, 800]);
    
    diff_data = [];
    diff_labels = {};
    
    for i = 1:n_levels
        level = intensity_levels{i};
        pre_vals = values_pre.(level);
        post_vals = values_post.(level);
        pre_vals = pre_vals(~isnan(pre_vals));
        post_vals = post_vals(~isnan(post_vals));
        
        % Compute difference (match lengths)
        min_len = min(length(pre_vals), length(post_vals));
        diff_vals = post_vals(1:min_len) - pre_vals(1:min_len);
        diff_data = [diff_data; diff_vals];
        diff_labels = [diff_labels; repmat({level}, length(diff_vals), 1)];
    end
    
    subplot(2, 3, 1);
    boxplot(diff_data, diff_labels);
    ylabel(sprintf('Difference (%s)', measure_name));
    title('Difference (postTMS - preTMS) across intensity levels');
    grid on;
    
    subplot(2, 3, 2);
    hold on;
    for i = 1:n_levels
        level = intensity_levels{i};
        pre_vals = values_pre.(level);
        post_vals = values_post.(level);
        pre_vals = pre_vals(~isnan(pre_vals));
        post_vals = post_vals(~isnan(post_vals));
        min_len = min(length(pre_vals), length(post_vals));
        diff_vals = post_vals(1:min_len) - pre_vals(1:min_len);
        histogram(diff_vals, 'FaceColor', colors(i,:), 'FaceAlpha', 0.6, ...
            'DisplayName', level);
    end
    xlabel(sprintf('Difference (%s)', measure_name));
    ylabel('Frequency');
    title('Histogram: Difference across intensity levels');
    legend('show');
    grid on;
    hold off;
    
    subplot(2, 3, 3);
    scatter_plot_multiple_groups_diff(values_pre, values_post, intensity_levels, measure_name);
    
    % Save second figure (differences)
    try
        fname2 = sprintf('%s_%s_differences.png', measure_type, measure_name);
        saveas(fig2, fullfile(outDir, fname2));
    catch
    end
    
    % Statistical tests
    fprintf('\nStatistical Tests:\n');
    fprintf('-------------------\n');
    
    % PreTMS vs postTMS within each level
    for i = 1:n_levels
        level = intensity_levels{i};
        pre_vals = values_pre.(level);
        post_vals = values_post.(level);
        pre_vals = pre_vals(~isnan(pre_vals));
        post_vals = post_vals(~isnan(post_vals));
        
        if length(pre_vals) > 0 && length(post_vals) > 0
            [h, p] = ttest2(pre_vals, post_vals);
            fprintf('%s: preTMS vs postTMS - p=%.4f (significant: %s)\n', ...
                level, p, char(string(h)));
        end
    end
    
    % PostTMS across levels
    if n_levels >= 2
        post_groups = cell(n_levels, 1);
        for i = 1:n_levels
            level = intensity_levels{i};
            post_vals = values_post.(level);
            post_groups{i} = post_vals(~isnan(post_vals));
        end
        [p_anova, ~] = anova1_multiple_groups(post_groups, intensity_levels);
        fprintf('postTMS across intensity levels - p=%.4f\n', p_anova);
    end
    
    % PreTMS across levels
    if n_levels >= 2
        pre_groups = cell(n_levels, 1);
        for i = 1:n_levels
            level = intensity_levels{i};
            pre_vals = values_pre.(level);
            pre_groups{i} = pre_vals(~isnan(pre_vals));
        end
        [p_anova, ~] = anova1_multiple_groups(pre_groups, intensity_levels);
        fprintf('preTMS across intensity levels - p=%.4f\n', p_anova);
    end
    
    % Difference across levels
    if n_levels >= 2
        diff_groups = cell(n_levels, 1);
        for i = 1:n_levels
            level = intensity_levels{i};
            pre_vals = values_pre.(level);
            post_vals = values_post.(level);
            pre_vals = pre_vals(~isnan(pre_vals));
            post_vals = post_vals(~isnan(post_vals));
            min_len = min(length(pre_vals), length(post_vals));
            diff_groups{i} = post_vals(1:min_len) - pre_vals(1:min_len);
        end
        [p_anova, ~] = anova1_multiple_groups(diff_groups, intensity_levels);
        fprintf('Difference (post-pre) across intensity levels - p=%.4f\n', p_anova);
    end
end

function scatter_plot_multiple_groups(values_struct, intensity_levels, condition, measure_name)
    hold on;
    colors = lines(length(intensity_levels));
    for i = 1:length(intensity_levels)
        level = intensity_levels{i};
        vals = values_struct.(level);
        vals = vals(~isnan(vals));
        scatter(1:length(vals), vals, 'filled', 'DisplayName', level, ...
            'MarkerFaceColor', colors(i,:));
    end
    xlabel('Sample');
    ylabel(measure_name);
    title(sprintf('Scatter: %s across intensity levels', condition));
    legend('show');
    grid on;
    hold off;
end

function scatter_plot_multiple_groups_diff(values_pre, values_post, intensity_levels, measure_name)
    hold on;
    colors = lines(length(intensity_levels));
    for i = 1:length(intensity_levels)
        level = intensity_levels{i};
        pre_vals = values_pre.(level);
        post_vals = values_post.(level);
        pre_vals = pre_vals(~isnan(pre_vals));
        post_vals = post_vals(~isnan(post_vals));
        min_len = min(length(pre_vals), length(post_vals));
        diff_vals = post_vals(1:min_len) - pre_vals(1:min_len);
        scatter(1:length(diff_vals), diff_vals, 'filled', 'DisplayName', level, ...
            'MarkerFaceColor', colors(i,:));
    end
    xlabel('Sample');
    ylabel(sprintf('Difference (%s)', measure_name));
    title('Scatter: Difference (post-pre) across intensity levels');
    legend('show');
    plot([0, max(cellfun(@length, struct2cell(values_pre)))], [0, 0], 'k--', 'LineWidth', 1);
    grid on;
    hold off;
end

function [p, table] = anova1_multiple_groups(groups, group_names)
%% Simple one-way ANOVA for multiple groups
    all_data = [];
    group_labels = {};
    
    for i = 1:length(groups)
        all_data = [all_data; groups{i}(:)];
        group_labels = [group_labels; repmat({group_names{i}}, length(groups{i}), 1)];
    end
    
    [p, table] = anova1(all_data, group_labels, 'off');
end
