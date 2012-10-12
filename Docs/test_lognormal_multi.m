function test_lognormal_multi()

    % Time params
    t_end = 60;
    dt = 0.1;
    t = 0:dt:t_end;

    % Geometry
    num_cells = 10;
    
    % Volume of water through
    flux_rate = 1;
    
    % Statistical params
    mu = 1;
    sigma = 0.5;
    
    % Initialize output
    out = zeros(num_cells, numel(t));
    
    % Modeling
    out(1, :) = flux_rate * dt * log_normal_pdf(t, mu, sigma, dt);
    for t_idx = 1:numel(t)
        for cell_idx = 2:num_cells
            out(cell_idx, t_idx:end) = out(cell_idx, t_idx:end) + out(cell_idx-1, t_idx) * dt * log_normal_pdf(t(t_idx:end) - t(t_idx), mu, sigma, dt);
        end
    end

    % Stochastic
    n_sim = 1e+4;
    out_st = zeros(num_cells, n_sim);
    
    rng_seed = 1;
    try
        rng(rng_seed);
    catch exception
        rand('seed', rng_seed);
    end
    
    out_st(1, :) = lognrnd(mu, sigma, 1, n_sim);
    for cell_idx = 2:num_cells
        out_st(cell_idx, :) = out_st(cell_idx-1, :) + lognrnd(mu, sigma, 1, n_sim);
    end

    % Plotting
    close all;
    num_hist_bars = 50;
    [n, xout] = hist(out_st(num_cells, :), num_hist_bars);
    factor = n_sim * (xout(2) - xout(1)) / (flux_rate * dt);
    bar(xout, n / factor, 1, 'r');
    
    hold on;
    plot(t, out, 'LineWidth', 2);
    legend('Stochastic', '1st', '2nd');    
    hold off;
    
    ev = mean(out_st, 2);
    stdev = std(out_st, 0, 2);
    figure(2);
%     plot(1:num_cells, ev, 'g');
%     hold on;
    plot(1:num_cells, (sqrt(stdev)), 'r');
%     hold off;
    
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