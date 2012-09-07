function execute_richards_1d_get_params

    EPSILON = 1e-8;

    close all;
    domain_name = 'channel';
%     domain_name = 'matrix';
    loop_type = 3;
    loop_type_names = {'ksat', 'wt', 'length'};

    data_wt = load(sprintf('../Common/data_%s_loop_%s_domain.mat', loop_type_names{loop_type}, domain_name));

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
        [mu_opt(i), sigma_opt(i), ratio_opt(i)] = pick_params(t, data_wt.out_flux(i, :));

        mu = mu_opt(i);
        sigma = sigma_opt(i);
        ratio = ratio_opt(i);

        out_flux_fickian = ratio * out_flux_lognrnd_pdf(t, mu, sigma);
        out_flux_fickian(1) = 0;

        abs_diff(i, :) = abs(data_wt.out_flux(i, :) - out_flux_fickian);
        
        disp(i / numel(var_vector) * 100);
        
        close all;
        plot(t, data_wt.out_flux(i, :));
        hold on;
        plot(t, out_flux_fickian, 'g');
        hold off;
    end

    params_opt = struct ();
    params_opt.saturation_effective_avg = data_wt.saturation_effective_avg;
    params_opt.k_sat_vector = data_wt.k_sat_vector;
    params_opt.van_genuchten_params = data_wt.van_genuchten_params;
    params_opt.mu = mu_opt;
    params_opt.sigma = sigma_opt;
    params_opt.ratio = ratio_opt;
    
    save(sprintf('../Common/opt_params_%s_%s_domain.mat', loop_type_names{loop_type}, domain_name), '-struct', 'params_opt');

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

    function [mux, sigmax, ratiox] = pick_params(t, out_flx, ratioi)
        mn = sum(out_flx .* t) / sum(out_flx);
        varn = sum(out_flx .* (t - mn).^2) / sum(out_flx);
        params = get_mu_sigma([mn varn]);
        mux = params(1);
        sigmax = params(2);
        out_flux_fickian = out_flux_lognrnd_pdf(t, mux, sigmax);
        if nargin < 3
            [~, ratiox] = square_diff(out_flx, out_flux_fickian);
        else
            ratiox = ratioi;
        end
    end

    function params = get_mu_sigma(moments)
        m = moments(1);
        v = moments(2);
        sigma = sqrt(log(v / (m * m) + 1));
        mu = log(m * exp(- sigma * sigma / 2));
        params = [mu sigma];
    end

    function res = out_flux_lognrnd_pdf(t, muv, sigmav)
        res = -1 ./ (sqrt(2 * pi) .* sigmav .* t) .* exp(-(log(t) - muv).^2 ./ (2 .* sigmav.^2));
        res(1) = 0;
    end

end