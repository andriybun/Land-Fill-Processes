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
    flux_rate = 10;
    
    % Statistical params
    log_normal_params = struct();
    log_normal_params.mu = [4; 1.5];
    log_normal_params.sigma = [0.6; 0.6];
    
    %% Parameters of domains
    domain_params = struct();
    domain_params.matrix_domain_idx = 1;
    domain_params.channel_domain_idx = 2;
    entry_ratio_matrix = 0.5;
    domain_params.entry_ratio = [entry_ratio_matrix; 1 - entry_ratio_matrix];
    
    tic;
    
    %% Baseline flux (log-normally distributed flux over all distance)
    t_idx = 1;
    out_flux_no_exchange = hydro(flux_rate, t_idx, geometry, domain_params, t_params, log_normal_params);
    out_flux_exchange = hydro_exchange(flux_rate, t_idx, geometry, domain_params, t_params, log_normal_params);
    
    t_idx = 300;
    out_flux_no_exchange = out_flux_no_exchange + hydro(flux_rate, t_idx, geometry, domain_params, t_params, log_normal_params);
    out_flux_exchange = out_flux_exchange + hydro_exchange(flux_rate, t_idx, geometry, domain_params, t_params, log_normal_params);
    
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

    return

    function out_flux = hydro(flux_rate, t_idx, geometry, domain_params, t_params, log_normal_params)
        %% Parsing input params
        t = t_params.t(t_idx:end) - t_params.t(t_idx);
        dt = t_params.dt;
        num_domains = geometry.num_domains;
        entry_ratio = domain_params.entry_ratio;
        mu = log_normal_params.mu;
        sigma = log_normal_params.sigma;
        
        %% Simulation
        out_flux = zeros(num_domains, numel(t_params.t));
        out_flux(1, t_idx:end) = entry_ratio(1) * flux_rate * dt * log_normal_pdf(t, mu(1), sigma(1), dt);
        out_flux(2, t_idx:end) = entry_ratio(2) * flux_rate * dt * log_normal_pdf(t, mu(2), sigma(2), dt);
    end
    
    function out_flux = hydro_exchange(flux_rate, t_idx, geometry, domain_params, t_params, log_normal_params)
        %% Parsing input params
        t = t_params.t(t_idx:end) - t_params.t(t_idx);
        dt = t_params.dt;
        num_cells = geometry.num_cells;
        num_domains = geometry.num_domains;
        matrix_domain_idx = domain_params.matrix_domain_idx;
        channel_domain_idx = domain_params.channel_domain_idx;
        entry_ratio = domain_params.entry_ratio;
        mu = log_normal_params.mu;
        sigma = log_normal_params.sigma;
    
        %% Calculating log-normal parameters for intermediate points:
        % Picked function
        pow = 1.7;
        steps = power(linspace(0, 1, num_cells + 1), pow);
        % Average flux time at the bottom (baseline)
        mean_orig = exp(mu + sigma .* sigma / 2);
        % Values of average flux time for intermediate points
        mean_steps = mean_orig * steps(2:end);
        % Mu for corresponding distributions of out-flux at intermediate points
        mu_steps = log(mean_steps) - repmat(sigma .* sigma / 2, [1, num_cells]);
        % Estimated mu's for response distributions from intermediate points to
        % the bottom
        mu_steps_out = log(repmat(mean_orig, [1, num_cells - 1]) - mean_steps(:, 1:end-1)) - repmat(sigma .* sigma / 2, [1, num_cells - 1]);
        
        %% Exchange rates function
        depth = 1:geometry.num_cells;
        exchange_rate_steps = zeros(geometry.num_domains, geometry.num_cells);
        exchange_rate_steps(matrix_domain_idx, :) = depth .^ 1.5;
        exchange_rate_steps(channel_domain_idx, 1) = 1;
        % Normalize:
        sm = repmat(sum(exchange_rate_steps, 2), [1, geometry.num_cells]);
        exchange_rate_steps = exchange_rate_steps ./ sm;
        
        %% Simulation
        out_flux = zeros(num_domains, numel(t_params.t));
        in_flux = flux_rate * repmat(exchange_rate_steps(:, 1) .* entry_ratio, [1, numel(t)]);
        out_flux(channel_domain_idx, t_idx:end) = out_flux(channel_domain_idx, t_idx:end) + ...
            in_flux(channel_domain_idx, :) .* dt .* log_normal_pdf(t, mu_steps(channel_domain_idx, num_cells), sigma(channel_domain_idx), dt);
        for cell_idx = 1:num_cells
            out_steps_channel_out = zeros(1, numel(t));
            in_flux = flux_rate * repmat(exchange_rate_steps(:, cell_idx) .* entry_ratio, [1, numel(t)]);
            out_flux_cell = in_flux(matrix_domain_idx, :) .* dt .* ...
                log_normal_pdf(t, mu_steps(matrix_domain_idx, cell_idx), sigma(matrix_domain_idx), dt);
            if cell_idx ~= num_cells
                for Xt_idx = 1:numel(t)
                    out_steps_dt_rep = repmat(out_flux_cell(Xt_idx) * dt, [1, numel(t) - Xt_idx + 1]);
                    out_steps_channel_out(Xt_idx:end) = out_steps_channel_out(Xt_idx:end) + ...
                        out_steps_dt_rep .* log_normal_pdf(t(Xt_idx:end) - t(Xt_idx), mu_steps_out(channel_domain_idx, cell_idx), sigma(channel_domain_idx), dt);
                end
            else
                out_flux(matrix_domain_idx, t_idx:end) = out_flux_cell;
            end
            out_flux(channel_domain_idx, t_idx:end) = out_flux(channel_domain_idx, t_idx:end) + out_steps_channel_out;
        end
    end
    
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