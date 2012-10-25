function inv_lognormal_multi()
% Log-normal transfer in 1D, one domain

    % Time params
    t_params = struct();
    t_end = 200;
    t_params.dt = 0.2;
    t_params.t = 0:t_params.dt:t_end;
    
    % Geometry
    geometry.num_cells = 2;
    num_cells =  geometry.num_cells;
    
    % Volume of water through
    flux_rate = 1;
    
    % Statistical params
    log_normal_params = struct();
    log_normal_params.mu = 2;
    log_normal_params.sigma = 1;
    
%     mu = log_normal_params.mu * ones(num_cells, 1);
%     sigma = log_normal_params.sigma * ones(num_cells, 1);
    
    % Initialize output
    out = zeros(num_cells, numel(t_params.t));
    
    % Modeling
    out(1, :) = hydro_1d(flux_rate, 1, geometry, t_params, log_normal_params);
    
    for t_idx = 1:numel(t_params.t)
        for cell_idx = 2:num_cells
            out(cell_idx, :) = out(cell_idx, :) + ...
                hydro_1d(out(cell_idx - 1, t_idx), t_idx, geometry, t_params, log_normal_params);
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
    
    out_st(1, :) = lognrnd(log_normal_params.mu, log_normal_params.sigma, 1, n_sim);
    out_st_mean(1) = mean(out_st(1, :));
    out_st_var(1) = var(out_st(1, :));
    for cell_idx = 2:num_cells
        tmp = lognrnd(log_normal_params.mu, log_normal_params.sigma, 1, n_sim);
        out_st(cell_idx, :) = out_st(cell_idx - 1, :) + tmp;
        out_st_mean(cell_idx) = mean(out_st(cell_idx, :));
        out_st_var(cell_idx) = var(out_st(cell_idx, :));
    end
    
    %% Plotting
    num_hist_bars = 50;
    [n, xout] = hist(out_st(num_cells, :), num_hist_bars);
    factor = n_sim * (xout(2) - xout(1)) / (flux_rate * t_params.dt);
    bar(xout, n / factor, 1, 'r');
    
    hold on;
    plot(t_params.t, out, 'LineWidth', 2);
    legend('Stochastic', '1st cell', '2nd cell');    
    hold off;
    
end