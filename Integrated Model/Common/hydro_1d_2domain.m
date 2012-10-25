function out_flux = hydro_1d_2domain(flux_rate, t_idx, geometry, domain_params, t_params, log_normal_params)
    %% Parsing input params
    t = t_params.t(t_idx:end) - t_params.t(t_idx);
    t = [t, t(end) + t_params.dt];
    dt = t_params.dt;
    num_domains = geometry.num_domains;
    entry_ratio = domain_params.entry_ratio;
    mu = log_normal_params.mu;
    sigma = log_normal_params.sigma;

    %% Simulation
    out_flux_cum = zeros(num_domains, numel(t_params.t) + 1);
    out_flux_cum(1, t_idx:end) = entry_ratio(1) * flux_rate * dt * log_normal_cdf(t, mu(1), sigma(1));
    out_flux_cum(2, t_idx:end) = entry_ratio(2) * flux_rate * dt * log_normal_cdf(t, mu(2), sigma(2));
    out_flux = out_flux_cum(:, 2:end) - out_flux_cum(:, 1:end-1);
end
