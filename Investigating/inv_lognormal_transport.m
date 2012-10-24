function test_lognormal_transport
    mu_ref = -3;
    z = 10;
    mu = 3.5075 * log(z) + mu_ref - 1;
    sigma = 0.6075;
    
    dt = 1e-1;
    t = 0:dt:500;
    out_flux = zeros(numel(t), numel(mu));
    
    for t_idx = 1:numel(t)
        out_flux(t_idx, :) = out_flux_lognrnd_pdf(t(t_idx), mu, sigma);
    end
    
    n_sim = 5000;
    out_flux_st = zeros(n_sim, numel(mu));
    for level_idx = 1:numel(mu)
        out_flux_st(:, level_idx) = lognrnd(mu(level_idx), sigma, n_sim, 1);
    end
    
%     close all
    
%     figure(1);
%     num_hist_bars = 50;
%     add_hist(out_flux_st(:, numel(mu)), num_hist_bars, 'r');
%     add_hist(out_flux_st(:, 1), num_hist_bars, 'g');
%     
%     sample = out_flux_st(:, 2) - out_flux_st(:, 1);
%     add_hist(sample, num_hist_bars, 'yellow');
    
    figure(3);
    hold on;
    plot(t, out_flux, 'LineWidth', 2, 'Color', 'k');
    hold off;
    
    min_val = 0.1;
    num_steps = 6;
%     steps = roundn(exp(linspace(log(min_val), log(1 + min_val), num_steps + 1)), -6) - min_val;

    pow = 1.7;
    steps = power(linspace(0, 1, num_steps + 1), pow);
    mean_orig = exp(mu + sigma * sigma / 2);
    mu_steps = log(mean_orig * steps(2:end)) - sigma * sigma / 2;
    
    out_flux_steps = zeros(numel(t), num_steps);
    for t_idx = 1:numel(t)
        out_flux_steps(t_idx, :) = out_flux_lognrnd_pdf(t(t_idx), mu_steps, sigma);
    end

    hold on;
    plot(t, out_flux_steps(:, 3:end));
    hold off;
%     [mu_calc, sigma_calc] = get_mu_sigma(mean(sample), var(sample));
%     mu_calc
%     sigma_calc
    
    return
    
    function res = out_flux_lognrnd_pdf(t, muv, sigmav)
        res = 1 ./ (sqrt(2 * pi) .* sigmav .* t) .* exp(-(log(t) - muv).^2 ./ (2 .* sigmav.^2));
        if t == 0
            res(:) = 0;
        end
    end

    function [mux, sigmax] = get_mu_sigma(evn, varn)
        sigmax = sqrt(log(varn ./ (evn .* evn) + 1));
        mux = log(evn) - sigmax .* sigmax ./ 2;
    end

    function add_hist(sample, num_hist_bars, color)
        hold on;
        num_sim = numel(sample);
        [n, xout] = hist(sample, num_hist_bars);
        factor = num_sim * (xout(2) - xout(1));
        bar(xout, n / factor, 1, color);
        hold off;
    end
end