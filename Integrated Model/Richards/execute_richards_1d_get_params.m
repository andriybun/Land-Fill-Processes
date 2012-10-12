function execute_richards_1d_get_params

    EPSILON = 1e-8;

    close all;
    
    domain_name = 'matrix';
%     domain_name = 'channel';
    
    loop_type = 3;
    loop_type_names = {'ksat', 'wt', 'length'};

    data_wt = load(sprintf('../Common/data_%s_loop_%s_domain.mat', loop_type_names{loop_type}, domain_name));
%     data_wt = load(sprintf('../Common/try_data_%s_domain.mat', domain_name));

    %% %%%%%%
    %%%
    %%%  Parameters analysis
    %%%
    %%%%%%%%%
    
    switch loop_type
        case 1          % saturated conductivity
            var_vector = data_wt.k_sat_vector;
        case 2          % water table
            var_vector = data_wt.water_table_elevation_vector;
            data_wt.k_sat_vector = repmat(data_wt.van_genuchten_params.k_sat, size(var_vector));
        case 3          % vG parameters
            var_vector = data_wt.length_vector;
            data_wt.k_sat_vector = repmat(data_wt.van_genuchten_params.k_sat, size(var_vector));
    end
    
    mu_opt = zeros(size(var_vector));
    sigma_opt = zeros(size(var_vector));
    ratio_opt = zeros(size(var_vector));
    
    abs_diff = zeros(size(data_wt.out_flux));
    
    for i = 1:numel(var_vector)
        t = data_wt.t_range(i, :);
        ratio = -sum(data_wt.out_flux(i, 1:end-1) + data_wt.out_flux(i, 2:end)) / 2 * (t(2) - t(1)) / 0.99;
        out_flux_normalized = data_wt.out_flux(i, :) / ratio;
        [mu, sigma] = pick_params(t, out_flux_normalized);

        mu_opt(i) = mu;
        sigma_opt(i) = sigma;
        ratio_opt(i) = ratio;

        dif_prev = inf;
        
        while true
            sigma = sigma + 2e-3;
            out_flux_lognorm = out_flux_lognrnd_pdf(t, mu, sigma);
            out_flux_lognorm(1) = 0;
            dif = square_diff(out_flux_normalized, out_flux_lognorm);
            if dif > dif_prev
                break
            end
            dif_prev = dif;
        end

        abs_diff(i, :) = abs(out_flux_normalized - out_flux_lognorm);
        
        disp(i / numel(var_vector) * 100);
        
        close all;
        figure('OuterPosition', [400, 600, 600, 520]);
        subplot(1, 1, 1);
        subplot('Position', [0.1 0.25 0.85 0.65]);

        plot(t, -out_flux_normalized, 'b');
        annotation('textbox', [0.1, 0.05, 0.85, 0.055], 'String', ...
            sprintf('S_e_f_f = %3.3f, K_s_a_t = %3.3f, theta_r = %3.3f, theta_s = %3.3f, alpha = %3.2f, lambda = %3.2f', ...
            data_wt.saturation_effective_avg(i),  data_wt.k_sat_vector(i), data_wt.van_genuchten_params.theta_r, ...
            data_wt.van_genuchten_params.theta_s, data_wt.van_genuchten_params.alpha, data_wt.van_genuchten_params.lambda));
        xlabel('Time');
        ylabel('Out flux');
        hold on;
        plot(t, -out_flux_lognorm, 'g');
        hold off;
        legend('Richard''s equation', sprintf('Log-normal (mu = %3.3f, sigma = %3.3f)', mu, sigma));
    end

    params_opt = struct ();
    params_opt.saturation_effective_avg = data_wt.saturation_effective_avg;
    params_opt.k_sat_vector = data_wt.k_sat_vector;
    params_opt.van_genuchten_params = data_wt.van_genuchten_params;
    params_opt.mu = mu_opt;
    params_opt.sigma = sigma_opt;
    params_opt.ratio = ratio_opt;
    
    save(sprintf('../Common/opt_params_%s_%s_domain.mat', loop_type_names{loop_type}, domain_name), '-struct', 'params_opt');
%     save(sprintf('../Common/try_opt_params_%s_domain.mat', domain_name), '-struct', 'params_opt');

    return
    
    function [res, ratio] = square_diff(v1, v2, ratiox)
        if nargin < 3
            mv1 = min(v1);
            mv2 = min(v2);
            ratio = mv1 / mv2;
        else
            ratio = ratiox;
        end
        res = sum((v1 - v2 .* ratio) .* (v1- v2 .* ratio));
    end

    function [mux, sigmax] = pick_params(t, out_flx)
        out_flx = -out_flx;
        dt = t(2:end) - t(1:end-1);
        mn = sum(t(2:end) .* (out_flx(1:end-1) + out_flx(2:end)) ./ 2 .* dt);
        varn = sum((t(2:end) - mn).^2 .* (out_flx(1:end-1) + out_flx(2:end)) ./ 2 .* dt);
        params = get_mu_sigma([mn varn]);
        mux = params(1);
        sigmax = params(2);
    end

    function params = get_mu_sigma(moments)
        m = moments(1);
        v = moments(2);
        sigmax = sqrt(log(v / (m * m) + 1));
        mux = log(m) - sigmax * sigmax / 2;
        params = [mux sigmax];
    end

    function res = out_flux_lognrnd_pdf(t, muv, sigmav)
        res = -1 ./ (sqrt(2 * pi) .* sigmav .* t) .* exp(-(log(t) - muv).^2 ./ (2 .* sigmav.^2));
        res(1) = 0;
    end

end