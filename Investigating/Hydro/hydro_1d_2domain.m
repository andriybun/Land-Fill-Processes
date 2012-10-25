function out_flux = hydro_1d_2domain(flux_rate, t_idx, geometry, domain_params, t_params, log_normal_params)
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
