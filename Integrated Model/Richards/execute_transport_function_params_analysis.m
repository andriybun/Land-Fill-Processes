function execute_transport_function_params_analysis()
    
    addpath('../Common/');

    domain_name = 'channel';
    domain_name = 'matrix';
    loop_type = 3;
    loop_type_names = {'ksat', 'wt', 'length'};

    data_wt = load(sprintf('../Common/data_%s_loop_%s_domain.mat', loop_type_names{loop_type}, domain_name));
    
%     opt_params_hc = load('opt_params_hc.mat');
    opt_params_wt = load(sprintf('../Common/opt_params_%s_%s_domain.mat', loop_type_names{loop_type}, domain_name));
    param_obj = log_normal_params(sprintf('../Common/opt_params_%s_%s_domain.mat', loop_type_names{2}, domain_name));
    
    close all;
    
    %% Simple analysis
    
    switch loop_type
        case 1
            var_vector = opt_params_wt.k_sat_vector;
            var_name = 'Hydraulic conductivity';
            const_val = opt_params_wt.saturation_effective_avg(1);
            const_name = 'S_e_f_f';
        case 2
            var_vector = opt_params_wt.saturation_effective_avg;
            var_name = 'Effective saturation';
            const_val = opt_params_wt.van_genuchten_params.k_sat;
            const_name = 'k_s_a_t';
        case 3
            var_vector = data_wt.length_vector;
            var_name = 'length';
            const_val = opt_params_wt.saturation_effective_avg(1);
            const_name = 'k_s_a_t = 1, S_e_f_f';
    end
    var_vector = roundn(var_vector, -8);
    
    %%
    plot(var_vector, opt_params_wt.mu);
    hold on;
    plot(var_vector, log(var_vector) * 2.5 - opt_params_wt.sigma .^ 2 ./ 2 + 2, 'g');
    hold off;
    %%
    
    figure('OuterPosition', [100, 100, 600, 400]);
    subplot(2, 2, 1);
    subplot('Position', [0.1 0.3 0.35 0.65]);
    plot(var_vector, opt_params_wt.mu);
%     if loop_type == 3
%         hold on;
%         plot(var_vector, opt_params_wt.saturation_effective_avg - 2, 'r');
%         hold off
%     end
    xlabel(var_name);
    ylabel('Mu');
    subplot(2, 2, 2);
    subplot('Position', [0.6 0.3 0.35 0.65]);
    plot(var_vector, opt_params_wt.sigma);
    xlabel(var_name);
    ylabel('Sigma');
    annotation('textbox', [0.1, 0.05, 0.85, 0.1], 'String', sprintf('%s = %3.3f, theta_r = %3.3f, theta_s = %3.3f, alpha = %3.2f, lambda = %3.2f', ...
        const_name, const_val, opt_params_wt.van_genuchten_params.theta_r, opt_params_wt.van_genuchten_params.theta_s, ...
        opt_params_wt.van_genuchten_params.alpha, opt_params_wt.van_genuchten_params.lambda));
    
%     %% Expression for mu
%     figure(2);
%     plot(var_vector, opt_params_wt.mu);
%     hold on;
%     v = -1.85;
%     plot(var_vector, log(1 ./ var_vector) + v, 'r');
%     hold off;
%     legend('mu', sprintf('log(1 / k_s_a_t) + %1.2f', v));
    
%     %% Water table loop
%     for i = 1:numel(data_wt.saturation_effective_avg)
%         t = data_wt.t_range(i, :);
%         [mu, sigm] = param_obj.get_params(data_wt.k_sat_vector(i), data_wt.saturation_effective_avg(i));
%         out_flux_approx = opt_params_wt.ratio(i) * ...
%             out_flux_lognrnd_pdf(t, mu, sigm);
%         close all;
%         hold on;
%         plot(t, data_wt.out_flux(i, :));
%         plot(t, out_flux_approx, 'g');
%         hold off;
%     end
    
    return

    function res = out_flux_lognrnd_pdf(t, muv, sigmav)
        res = -1 ./ (sqrt(2 * pi) .* sigmav .* t) .* exp(-(log(t) - muv).^2 ./ (2 .* sigmav.^2));
        res(1) = 0;
    end

end