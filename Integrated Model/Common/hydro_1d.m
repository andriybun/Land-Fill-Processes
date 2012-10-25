function out_flux = hydro_1d(flux_rate, t_idx, geometry, t_params, log_normal_params)
    %% Parsing input params
    t = t_params.t(t_idx:end) - t_params.t(t_idx);
    t = [t, t(end) + t_params.dt];
    dt = t_params.dt;
    mu = log_normal_params.mu;
    sigma = log_normal_params.sigma;

    %% Simulation
    out_flux_cum = zeros(1, numel(t_params.t) + 1);
    out_flux_cum(t_idx:end) = flux_rate * dt * log_normal_cdf(t, mu, sigma);
    out_flux = out_flux_cum(2:end) - out_flux_cum(1:end-1);
end
