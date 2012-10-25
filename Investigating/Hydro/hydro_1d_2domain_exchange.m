function out_flux = hydro_1d_2domain_exchange(flux_rate, t_idx, geometry, domain_params, t_params, log_normal_params)
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
    % Picked function that defines distribution of mean average flow over cells
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
    % Function that defines how much flux swithces from matrix to channel
    % domain after each cell
    exchange_rate_steps(matrix_domain_idx, :) = 1; %depth .^ 1.5;
    % All water that enters channel domain flows out directly
    exchange_rate_steps(channel_domain_idx, 1) = 1;
    % Normalize:
    sm = repmat(sum(exchange_rate_steps, 2), [1, geometry.num_cells]);
    exchange_rate_steps = exchange_rate_steps ./ sm;

    %% Simulation
    out_flux = zeros(num_domains, numel(t_params.t));
    in_flux = flux_rate * repmat(exchange_rate_steps(:, 1) .* entry_ratio, [1, numel(t)]);
    % Output flux of water that initially entered channel domain
    out_flux(channel_domain_idx, t_idx:end) = out_flux(channel_domain_idx, t_idx:end) + ...
        in_flux(channel_domain_idx, :) .* dt .* log_normal_pdf(t, mu_steps(channel_domain_idx, num_cells), sigma(channel_domain_idx), dt);
    %
    for cell_idx = 1:num_cells
        out_steps_channel_out = zeros(1, numel(t));
        in_flux = flux_rate * repmat(exchange_rate_steps(:, cell_idx) .* entry_ratio, [1, numel(t)]);
        % Amount of water that is transferred from matrix to channel domain
        out_flux_cell = in_flux(matrix_domain_idx, :) .* dt .* ...
            log_normal_pdf(t, mu_steps(matrix_domain_idx, cell_idx), sigma(matrix_domain_idx), dt);
        if cell_idx ~= num_cells
            for Xt_idx = 1:numel(t)
                out_steps_dt_rep = repmat(out_flux_cell(Xt_idx) * dt, [1, numel(t) - Xt_idx + 1]);
                % Calculate outflux of water that has been transferred
                % from matrix to channel domain
                out_steps_channel_out(Xt_idx:end) = out_steps_channel_out(Xt_idx:end) + ...
                    out_steps_dt_rep .* log_normal_pdf(t(Xt_idx:end) - t(Xt_idx), mu_steps_out(channel_domain_idx, cell_idx), sigma(channel_domain_idx), dt);
            end
        else
            % Water that remains in matrix domain until bottom flows out
            out_flux(matrix_domain_idx, t_idx:end) = out_flux_cell;
        end
        out_flux(channel_domain_idx, t_idx:end) = out_flux(channel_domain_idx, t_idx:end) + out_steps_channel_out;
    end
end