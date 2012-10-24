function test_lognormal_multi()

    % Time params
    t_end = 100;
    dt = 0.1;
    t = 0:dt:t_end;

    % Geometry
    num_cells = 2;
    
    % Volume of water through
    flux_rate = 1;
    
    % Statistical params
    mu = 2;
    sigma = 1;
    
    mu = repmat(mu, [num_cells 1]);
    sigma = repmat(sigma, [num_cells 1]);
    
    % Initialize output
    out = zeros(num_cells, numel(t));
    
    % Modeling
    out(1, :) = flux_rate * dt * log_normal_pdf(t, mu(1), sigma(1), dt);
    for t_idx = 1:numel(t)
        for cell_idx = 2:num_cells
            out(cell_idx, t_idx:end) = out(cell_idx, t_idx:end) + out(cell_idx - 1, t_idx) * dt * ...
                log_normal_pdf(t(t_idx:end) - t(t_idx), mu(cell_idx), sigma(cell_idx), dt);
        end
    end

    % Stochastic
    n_sim = 1e+4;
    out_st = zeros(num_cells, n_sim);
    out_st_mean = zeros(num_cells, 1);
    out_st_var = zeros(num_cells, 1);
    
    rng_seed = 1;
    try
        rng(rng_seed);
    catch dummy
        rand('seed', rng_seed);
    end
    
    out_st(1, :) = lognrnd(mu(1), sigma(1), 1, n_sim);
    out_st_mean(1) = mean(out_st(1, :));
    out_st_var(1) = var(out_st(1, :));
    for cell_idx = 2:num_cells
        tmp = lognrnd(mu(cell_idx), sigma(cell_idx), 1, n_sim);
        out_st(cell_idx, :) = out_st(cell_idx - 1, :) + tmp;
        out_st_mean(cell_idx) = mean(out_st(cell_idx, :));
        out_st_var(cell_idx) = var(out_st(cell_idx, :));
    end
    
    %% Plotting
    num_hist_bars = 50;
    [n, xout] = hist(out_st(num_cells, :), num_hist_bars);
    factor = n_sim * (xout(2) - xout(1)) / (flux_rate * dt);
    bar(xout, n / factor, 1, 'r');
    
    hold on;
    plot(t, out, 'LineWidth', 2);
    legend('Stochastic', '1st cell', '2nd cell');    
    hold off;
    
%     [cumsum(out_st_var), var(out_st')']
%     [cumsum(out_st_mean), mean(out_st')']
    
    return

    function res = log_normal_pdf(t_in, mu_in, sigma_in, dt)
        t_in = repmat(t_in, [numel(mu_in), 1]);
        mu_in = repmat(mu_in, [1, size(t_in, 2)]);
        sigma_in = repmat(sigma_in, [1, size(t_in, 2)]);
        res = exp(-(log(t_in) - mu_in).^2 ./ (2 .* sigma_in .* sigma_in)) ./ (sqrt(2 * pi) .* sigma_in .* t_in);
        sigma_is_inf = isinf(sigma_in(:, 1));
        res(~sigma_is_inf, 1) = 0;
        res(sigma_is_inf, 1) = 1 / dt;
    end
end