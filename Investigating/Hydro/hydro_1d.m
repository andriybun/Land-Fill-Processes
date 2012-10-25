function out_flux = hydro_1d(flux_rate, t_idx, geometry, t_params, log_normal_params)
    %% Parsing input params
    t = t_params.t(t_idx:end) - t_params.t(t_idx);
    dt = t_params.dt;
    mu = log_normal_params.mu;
    sigma = log_normal_params.sigma;

    %% Simulation
    out_flux = zeros(1, numel(t_params.t));
    out_flux(t_idx:end) = flux_rate * dt * log_normal_pdf(t, mu, sigma, dt);
end
