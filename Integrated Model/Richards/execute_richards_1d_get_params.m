function execute_richards_1d_get_params

    EPSILON = 1e-8;

    close all;
    domain_name = 'matrix';
    data_wt = load(sprintf('../Common/data_wt_loop_%s_domain.mat', domain_name));
%     data_hc = load('data_hydro_conductivity_loop.mat');

	%% Compare breakthrough and effective saturation
%     alpha   = 0.9;
% 	  theta_s = 0.368;
%     theta_r = 0.102;
%     lambda  = 2.1;
%     n       = lambda + 1;
%     m       = lambda / n;
%     
%     h = -2.5:0.1:1;
%     hc = -h;
%     se = (1 + abs(alpha .* hc).^n).^(-m) .* (hc > 0) + 1 .* (hc <= 0);
%     se = 1 - (1 - se) .* 0.78 ./ 0.85;
%     
%     % Compute the effective saturation
%     theta = se .* (theta_s - theta_r) + theta_r;
% 
%     hold on;
%     plot(data_wt.water_table_elevation_vector, data_wt.saturation_effective_avg, 'b');
%     plot(h, se, 'r');
%     hold off;
    
    %% 
    
    t = data_wt.t_range;
    
    %% %%%%%%
    %%%
    %%%  Water table analysis
    %%%
    %%%%%%%%%
    
    z_opt = zeros(size(data_wt.water_table_elevation_vector));
    v_opt = zeros(size(data_wt.water_table_elevation_vector));
    d_opt = zeros(size(data_wt.water_table_elevation_vector));
    ratio_opt = zeros(size(data_wt.water_table_elevation_vector));
    
    abs_diff = zeros(size(data_wt.out_flux));
    
    for i = 1:numel(data_wt.water_table_elevation_vector)
        t = data_wt.t_range(i, :);
        [z_opt(i), v_opt(i), d_opt(i), ratio_opt(i)] = pick_params(t, data_wt.out_flux(i, :), data_wt.water_table_elevation_vector(i));

        z = z_opt(i);
        v = v_opt(i);
        d = d_opt(i);
        ratio = ratio_opt(i);

        out_flux_fickian = ratio * out_flux_lognrnd_pdf(t, z, v, d);
        out_flux_fickian(1) = 0;

        abs_diff(i, :) = abs(data_wt.out_flux(i, :) - out_flux_fickian);
        
        disp(i / numel(data_wt.water_table_elevation_vector) * 100);
        
        close all;
        plot(t, data_wt.out_flux(i, :));
        hold on;
        plot(t, out_flux_fickian, 'g');
        hold off;
    end

    
%     idx = ((abs(data_wt.out_flux) >= EPSILON) & (abs_diff >= 1e+4 * EPSILON));
%     rel_err = zeros(size(data_wt.out_flux));
%     rel_err(idx) = abs_diff(idx) ./ abs(data_wt.out_flux(idx)) * 100;
%     mesh(rel_err);
%     disp(max(max(rel_err)));

    params_opt_wt = struct ();
    params_opt_wt.saturation_effective_avg = data_wt.saturation_effective_avg;
    params_opt_wt.van_genuchten_params = data_wt.van_genuchten_params;
    params_opt_wt.z = z_opt;
    params_opt_wt.v = v_opt;
    params_opt_wt.d = d_opt;
    params_opt_wt.ratio = ratio_opt;
    
    save(sprintf('../Common/opt_params_wt_%s_domain.mat', domain_name), '-struct', 'params_opt_wt');

    %% %%%%%%
    %%%
    %%%  Hydraulic conductivity analysis
    %%%
    %%%%%%%%%
    
    t = data_hc.t_range;
    
    z_opt = zeros(size(data_hc.k_sat_vector));
    v_opt = zeros(size(data_hc.k_sat_vector));
    d_opt = zeros(size(data_hc.k_sat_vector));
    ratio_opt = zeros(size(data_hc.k_sat_vector));
    
    for i = 1:numel(data_hc.k_sat_vector)
        t = data_hc.t_range(i, :);
        [z_opt(i), v_opt(i), d_opt(i), ratio_opt(i)] = pick_params(t, data_hc.out_flux(i, :), data_hc.k_sat_vector(i));

        z = z_opt(i);
        v = v_opt(i);
        d = d_opt(i);
        ratio = ratio_opt(i);
        
        disp(i / numel(data_hc.k_sat_vector) * 100);
    end
    
    params_opt_hc = struct ();
    params_opt_hc.k_sat_vector = data_hc.k_sat_vector;
    params_opt_hc.z = z_opt;
    params_opt_hc.v = v_opt;
    params_opt_hc.d = d_opt;
    params_opt_hc.ratio = ratio_opt;
    
    save('opt_params_hc.mat', '-struct', 'params_opt_hc');

    copyfile('opt_params_wt.mat', '~/Dropbox/');
    copyfile('opt_params_hc.mat', '~/Dropbox/');
    
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

    function [zx, vx, dx, ratiox] = pick_params(t, out_flx, var, ratioi)
        mn = sum(out_flx .* t) / sum(out_flx);
        varn = sum(out_flx .* (t - mn).^2) / sum(out_flx);
        params = get_mu_sigma([mn varn]);
        zx = params(1); % 3.4 * log(-var) - 0.24; % -log(var) - 0.3915658;
        vx = params(2); % 0.627; % 
        dx = 0;
        out_flux_fickian = out_flux_lognrnd_pdf(t, zx, vx, dx);
        if nargin < 4
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
%     function moments = get_mu_sigma(params)
%         mu = params(1);
%         sigma = params(2);
%         moments = [
%             exp(mu + sigma * sigma / 2)
%             (exp(sigma * sigma) - 1) * exp(2 * mu + sigma * sigma);
%         ];
%     end

    function res = out_flux_fickian_pdf(t, zv, vv, dv)
        res = -zv .* exp(-(zv - vv .* t).^2 ./ (4 .* dv .* t)) ./ (2 .* sqrt(pi .* dv .* t.^3));
        res(1) = 0;
    end

    function res = out_flux_lognrnd_pdf(t, muv, sigmav, dummy)
        res = -1 ./ (sqrt(2 * pi) .* sigmav .* t) .* exp(-(log(t) - muv).^2 ./ (2 .* sigmav.^2));
        res(1) = 0;
    end
    
    function res = z_wt(wt)
        res = - 80 * wt;
    end

    function res = v_wt(wt)
        res = 7.19 * (-wt) ^ -2.49;
    end

    function res = d_wt(wt)
        res = 202.3 * (-wt) ^ -1.36;
    end

    function res = z_hc(k_satur)
        res = 15.8 * k_satur ^ -0.6;
    end

    function res = v_hc(k_satur)
        res = 15.46 * k_satur ^ 0.37;
    end

    function res = d_hc(k_satur)
        res = 80.565 * k_satur ^ -0.22;
    end
end