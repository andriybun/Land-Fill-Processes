function inv_lognormal_multi_2()

    %% Inputs
    
    % Time params
    t_params = struct();
    t_end = 200;
    t_params.dt = 0.2;
    t_params.t = 0:t_params.dt:t_end;

    % geometryetry
    geometry = struct();
    geometry.num_cells = 5;
    geometry.num_domains = 2;
    
    % Volume of water through
    flux_rate = 1;
    
    % Statistical params
    log_normal_params = struct();
    log_normal_params.mu = [4; 1.5];
    log_normal_params.sigma = [0.6; 0.6];
    
    %% Parameters of domains
    domain_params = struct();
    domain_params.num_domains = 2;
    domain_params.matrix_domain_idx = 1;
    domain_params.channel_domain_idx = 2;
    % Fraction of water that initially enters matrix domain
    entry_ratio_matrix = 0.9;
    domain_params.entry_ratio = [entry_ratio_matrix; 1 - entry_ratio_matrix];
    
    tic;
    
    %% Baseline flux (log-normally distributed flux over all distance)
    t_idx = 1;
    out_flux_no_exchange = hydro_1d_2domain(flux_rate, t_idx, geometry, domain_params, t_params, log_normal_params);
    out_flux_exchange = hydro_1d_2domain_exchange(flux_rate, t_idx, geometry, domain_params, t_params, log_normal_params);
    
    t_idx = 300;
    out_flux_no_exchange = out_flux_no_exchange + hydro_1d_2domain(flux_rate, t_idx, geometry, domain_params, t_params, log_normal_params);
    out_flux_exchange = out_flux_exchange + hydro_1d_2domain_exchange(flux_rate, t_idx, geometry, domain_params, t_params, log_normal_params);
    
    toc;
    
    %% Plotting
    close all;
    hold on;
    plot(t_params.t, squeeze(sum(out_flux_exchange, 1)), 'LineWidth', 2, 'Color', [0.9, 0, 0.9]);
    plot(t_params.t, sum(out_flux_no_exchange, 1), 'LineWidth', 2, 'Color', [0.9, 0.9, 0]);
    legend('simulated', 'baseline');
    hold off;

    %% Checking mass balance
    break_steps = sum(out_flux_exchange, 2)
    disp(sum(out_flux_no_exchange, 2))

end