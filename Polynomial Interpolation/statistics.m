function statistics(func, a, b, n)
    % Define node types and orderings
    node_types = {'uniform', 'chebyshev1', 'chebyshev2'};
    barycentric2_flags = {'uniform', 'cheby1', 'cheby2'};
    orderings = {'ascend', 'descend', 'leja'};
    order_colors = {'r', 'b', 'k'};  % Red, Blue, Black for orders

    % Initialize arrays to store results
    results = [];
    conditioning_results = [];

    % Define test points
    xx = linspace(a, b, 100);
    ff = func(xx);  

    % Loop over node types
    for i = 1:length(node_types)
        node_type = node_types{i};
        bary2_flag = barycentric2_flags{i};

        % Generate base nodes
        base_nodes = order_function(node_type, a, b, n);

        % Compute conditioning measures
        [gamma, f_bary1] = barycentric1(func, n, base_nodes);
        [~, ~, kappa_one] = bary1Eval(xx, f_bary1, base_nodes, gamma);
        lebesgue_const = max(kappa_one);

        [beta, f_bary2] = barycentric2(func, n, bary2_flag);
        ff_bary2 = bary2Eval(xx, f_bary2, base_nodes, beta, n);
        script_H_n = max(abs(ff_bary2 - ff));

        % Store conditioning results
        conditioning_results = [conditioning_results; {node_type, lebesgue_const, script_H_n}];

        % Loop over orderings
        for j = 1:length(orderings)
            ordering = orderings{j};
            nodes = order(base_nodes, ordering);

            % Compute interpolations and errors
            [dd, f_newton] = newtondd(func, nodes);
            ff_newton = horner(dd, nodes, xx);
            abs_error_newton = max(abs(ff_newton - ff));

            [gamma, f_bary1] = barycentric1(func, n, nodes);
            [ff_bary1, ~, ~] = bary1Eval(xx, f_bary1, nodes, gamma);
            abs_error_bary1 = max(abs(ff_bary1 - ff));

            [beta, f_bary2] = barycentric2(func, n, bary2_flag);
            ff_bary2 = bary2Eval(xx, f_bary2, nodes, beta, n);
            abs_error_bary2 = max(abs(ff_bary2 - ff));

            % Store results
            results = [results; {func2str(func), node_type, ordering, abs_error_bary1, abs_error_bary2, abs_error_newton}];
        end
    end

    % Convert conditioning results to a LaTeX table
    generate_latex_table(conditioning_results, func2str(func));

    % Plot histograms with colored bars
    plot_separate_histograms(results, func2str(func), node_types, orderings, order_colors);
end

% Function to determine which node function to call
function nodes = order_function(node_type, a, b, n)
    switch node_type
        case 'uniform'
            nodes = uniform(a, b, n);
        case 'chebyshev1'
            nodes = chebyshev1(a, b, n);
        case 'chebyshev2'
            nodes = chebyshev2(a, b, n);
        otherwise
            error('Invalid node type.');
    end
end
% Function to generate LaTeX table for conditioning metrics
function generate_latex_table(conditioning_results, func_name)
    fprintf('\\begin{table}[h]\n');
    fprintf('\\centering\n');
    fprintf('\\begin{tabular}{|c|c|c|}\n');
    fprintf('\\hline\n');
    fprintf('Node Type & Lebesgue Constant & H_n \\\\\n');
    fprintf('\\hline\n');
    
    for i = 1:size(conditioning_results, 1)
        fprintf('%s & %.4f & %.4f \\\\\n', ...
            conditioning_results{i, 1}, ...
            conditioning_results{i, 2}, ...
            conditioning_results{i, 3});
    end
    
    fprintf('\\hline\n');
    fprintf('\\end{tabular}\n');
    fprintf('\\caption{Conditioning Metrics for %s}\n', func_name);
    fprintf('\\label{tab:%s-conditioning}\n', func_name);
    fprintf('\\end{table}\n\n');
end

% Function to plot separate histograms with colored bars per ordering
function plot_separate_histograms(results, func_name, node_types, orderings, order_colors)
    % Convert results to table for easier indexing
    results_table = cell2table(results, 'VariableNames', ...
        {'Function', 'NodeType', 'Ordering', 'MaxAbsError_Bary1', 'MaxAbsError_Bary2', 'MaxAbsError_Newton'});

    % Extract data grouped by method, then by node type and ordering
    errors_bary1 = zeros(length(node_types), length(orderings));
    errors_bary2 = zeros(length(node_types), length(orderings));
    errors_newton = zeros(length(node_types), length(orderings));

    for i = 1:length(node_types)
        for j = 1:length(orderings)
            row = strcmp(results_table.NodeType, node_types{i}) & strcmp(results_table.Ordering, orderings{j});
            errors_bary1(i, j) = results_table.MaxAbsError_Bary1(row);
            errors_bary2(i, j) = results_table.MaxAbsError_Bary2(row);
            errors_newton(i, j) = results_table.MaxAbsError_Newton(row);
        end
    end

    % X-axis labels
    x_labels = node_types;

    % Histogram for Newton
    figure;
    hold on;
    for j = 1:length(orderings)
        bar((1:length(node_types)) + (j-2)*0.2, errors_newton(:, j), 0.2, order_colors{j}, 'DisplayName', orderings{j});
    end
    hold off;
    set(gca, 'XTick', 1:length(x_labels), 'XTickLabel', x_labels, 'XTickLabelRotation', 45);
    xlabel('Node Type');
    ylabel('Max Absolute Error');
    title(['Newton Interpolation Error for ', func_name]);
    legend show;
    grid on;

    % Histogram for Barycentric 1
    figure;
    hold on;
    for j = 1:length(orderings)
        bar((1:length(node_types)) + (j-2)*0.2, errors_bary1(:, j), 0.2, order_colors{j}, 'DisplayName', orderings{j});
    end
    hold off;
    set(gca, 'XTick', 1:length(x_labels), 'XTickLabel', x_labels, 'XTickLabelRotation', 45);
    xlabel('Node Type');
    ylabel('Max Absolute Error');
    title(['Barycentric 1 Interpolation Error for ', func_name]);
    legend show;
    grid on;

    % Histogram for Barycentric 2
    figure;
    hold on;
    for j = 1:length(orderings)
        bar((1:length(node_types)) + (j-2)*0.2, errors_bary2(:, j), 0.2, order_colors{j}, 'DisplayName', orderings{j});
    end
    hold off;
    set(gca, 'XTick', 1:length(x_labels), 'XTickLabel', x_labels, 'XTickLabelRotation', 45);
    xlabel('Node Type');
    ylabel('Max Absolute Error');
    title(['Barycentric 2 Interpolation Error for ', func_name]);
    legend show;
    grid on;
end
