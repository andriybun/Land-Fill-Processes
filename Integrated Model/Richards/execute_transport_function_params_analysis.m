function execute_transport_function_params_analysis()
    
%     data_hc = load('data_hydro_conductivity_loop.mat');
    domain_name = 'matrix';
    data_wt = load(sprintf('../Common/data_wt_loop_%s_domain.mat', domain_name));
    
%     opt_params_hc = load('opt_params_hc.mat');
    opt_params_wt = load(sprintf('../Common/opt_params_wt_%s_domain.mat', domain_name));
    param_obj = log_normal_params(sprintf('../Common/opt_params_wt_%s_domain.mat', domain_name));
    
    close all;
    
    %% Water table loop
    for i = 1:numel(data_wt.saturation_effective_avg)
        t = data_wt.t_range(i, :);
        [mu, sigm] = param_obj.get_params(1, data_wt.saturation_effective_avg(i));
        out_flux_approx = opt_params_wt.ratio(i) * ...
            out_flux_lognrnd_pdf(t, mu, sigm, 0);
        close all;
        hold on;
        plot(t, data_wt.out_flux(i, :));
        plot(t, out_flux_approx, 'g');
        hold off;
    end
    
%     %% General
%     opt_params = opt_params_wt;
%  
%     figure(2);
%     hold on;
%     plot(data_wt.water_table_elevation_vector(4:end), opt_params.v(4:end), '-g*');
%     plot(data_wt.water_table_elevation_vector(4:end), opt_params.z(4:end), '-b*');
% %     plot(data_wt.water_table_elevation_vector(4:end), 1e+5 * opt_params.ratio(4:end), '-k*');
% %     plot(data_wt.saturation_effective_avg(4:end), opt_params.v(4:end), '-g*');
% %     plot(data_wt.saturation_effective_avg(4:end), opt_params.z(4:end), '-b*');
% %     plot(data_wt.saturation_effective_avg(4:end), 1e+5 * opt_params.ratio(4:end), '-k*');
%     hold off;

    %% Hydraulic conductivity loop
    for i = 1:numel(data_hc.k_sat_vector)
        t = data_hc.t_range(i, :);
%         out_flux_fickian = opt_params_hc.ratio(i) * ...
%             out_flux_lognrnd_pdf(t, opt_params_hc.z(i), opt_params_hc.v(i), opt_params_hc.d(i));
        [mu, sigm] = param_obj.get_params(data_hc.k_sat_vector(i), 0.3255);
        out_flux_fickian = opt_params_hc.ratio(i) * ...
            out_flux_lognrnd_pdf(t, mu, sigm, 0);
        close all;
        hold on;
        plot(t, data_hc.out_flux(i, :));
        plot(t, out_flux_fickian, 'g');
        hold off;
    end

%     %% General
%     opt_params = opt_params_hc;
%  
%     figure(2);
%     hold on;
%     plot(data_hc.k_sat_vector, opt_params.v, '-g*');
%     plot(data_hc.k_sat_vector, opt_params.z, '-b*');
%     plot(data_hc.k_sat_vector, 1e+5 * opt_params.ratio, '-k*');
%     hold off;

    return
    
    function res = out_flux_fickian_pdf(t, zv, vv, dv)
        res = -zv .* exp(-(zv - vv .* t).^2 ./ (4 .* dv .* t)) ./ (2 .* sqrt(pi .* dv .* t.^3));
        res(1) = 0;
    end

    function res = out_flux_lognrnd_pdf(t, muv, sigmav, dummy)
        res = -1 ./ (sqrt(2 * pi) .* sigmav .* t) .* exp(-(log(t) - muv).^2 ./ (2 .* sigmav.^2));
        res(1) = 0;
    end

    function res = se(hc, vg_par)
        res = (1 + (vg_par.alpha .* hc).^vg_par.n).^(-vg_par.m) .* (hc > 0) + 1 .* (hc <= 0);
    end

    % Sigma of effective saturation
    function y = sigma(se)
        p1 = 4.4014;
        p2 = -9.8347;
        p3 = 8.0144;
        p4 = -2.6853;
        p5 = 0.20365;
        p6 = 0.66483;
        
        y = p1 * se^5 + p2 * se^4 + p3 * se^3 + p4 * se^2 + p5 * se + p6;
    end

end