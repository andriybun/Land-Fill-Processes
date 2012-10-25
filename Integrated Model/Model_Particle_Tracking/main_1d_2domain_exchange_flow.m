function main_1d_2domain_exchange_flow()
    clc;
    tic;
    addpath('../Common/');

%     % Time params
%     start_date = struct();
%     start_date.year = 2012;
%     start_date.month = 1;
%     start_date.day = 1;
%     t_params.max_days = 100;                                                           % number of simulation days
%     t_params.intervals_per_day = 24;
%     t_params.dt = 1 / t_params.intervals_per_day;                 % in days
%     t_params.num_intervals = t_params.max_days * t_params.intervals_per_day;    % in {time step}
%     t_params.days_elapsed = (0 : 1: (t_params.num_intervals-1)) / t_params.intervals_per_day;
%     num_intervals = t_params.num_intervals;
    
    % Geometry params
    geometry_params = struct();
    geometry_params.dz = 1;
    geometry_params.zn = 10;
    geometry_params.is_landfill_array = ones(geometry_params.zn, 1);
    geometry_params.column_height_array = geometry_params.dz * squeeze(sum(geometry_params.is_landfill_array, 1));
    geometry_params.num_domains = 2;

    % Precipitation data
    file_name = '../Data/precipitation_daily_KNMI_20110908.txt';
    [precipitation_intensity_time_vector, t_params, start_date] = read_precipitation_data_csv(file_name);
    num_intervals = t_params.num_intervals;
%     %% Input pulse:
%     precipitation_intensity_time_vector = zeros(size(t_params.days_elapsed));
%     precipitation_intensity_time_vector(:) = 0;
%     precipitation_intensity_time_vector(1) = 1e-4 * t_params.dt;
%     %% End input pulse

    %% TESTING (let all the water flow out to check mass balance):
    extra_days = 100;
    t_params.max_days = t_params.max_days + extra_days;
    t_params.num_intervals = t_params.max_days * t_params.intervals_per_day;                           % in {time step}
    num_intervals = t_params.num_intervals;
    t_params.t = 0:t_params.dt:t_params.max_days - t_params.dt;
    t_params.days_elapsed = (0 : 1: (t_params.num_intervals-1)) / t_params.intervals_per_day;
    precipitation_intensity_time_vector = cat(2, precipitation_intensity_time_vector, zeros(1, extra_days * t_params.intervals_per_day));
    %% END TESTING

    %% Parameters of domains
    domain_params = struct();
    domain_params.matrix_domain_idx = 1;
    domain_params.channel_domain_idx = 2;
    % Fraction of water that initially enters matrix domain
    entry_ratio_matrix = 0.5;
    domain_params.entry_ratio = [entry_ratio_matrix; 1 - entry_ratio_matrix];

    % Determine probability distribution parameters corresponding to defined inputs:
    lognrnd_param_definer(1) = log_normal_params('../Common/opt_params_wt_matrix_domain.mat');
    lognrnd_param_definer(2) = log_normal_params('../Common/opt_params_wt_channel_domain.mat');
    
    % Fluid velocity parameters (1st - matrix domain; 2nd - channel domain)
    % Note: vectors must be row vectors
    hydraulic_params            = aos_to_soa([lognrnd_param_definer(1).hydraulic_params
                                              lognrnd_param_definer(2).hydraulic_params]');
    hydraulic_params.k_sat_ref  = hydraulic_params.k_sat;
    hydraulic_params.k_sat      = [1, 1];
	hydraulic_params.d          = [1, 1];                                                        % diffusion_coefficient

    % initial SE
    properties_array = generate_biogeochemical_properties(geometry_params, hydraulic_params);
    
    % Result
    leachate_flux = zeros(geometry_params.num_domains, num_intervals);
    
    progr = floor(t_params.num_intervals / 10);
    
    for t_idx = 1:num_intervals
        [leachate_flux, properties_array] = transport_lognormal(leachate_flux, t_idx, ...
            properties_array, precipitation_intensity_time_vector(t_idx), geometry_params, ...
            domain_params, hydraulic_params, t_params, lognrnd_param_definer);
        
        % Display progress
        if (mod(t_idx, progr) == 0)
            disp (t_idx / t_params.num_intervals * 100);
        end
    end

    plot(t_params.days_elapsed, sum(leachate_flux, 3));
    plot(t_params.days_elapsed, squeeze(leachate_flux(1, :)), 'r');
    hold on;
    plot(t_params.days_elapsed, squeeze(leachate_flux(2, :)), 'g');
    plot(t_params.days_elapsed, squeeze(sum(leachate_flux, 1)), 'k');
    hold off;
%     hold on;
%     plot(t_params.days_elapsed, leachate_flux_matrix_domain, 'g');
%     plot(t_params.days_elapsed, leachate_flux_fast_domain, 'r');
%     hold off;
%     figure(2);
%     mesh(se);

    toc;

    return 
    
    function [leachate_intercell_array, properties_array] = transport_lognormal(leachate_intercell_array, t_idx, properties_array, scale, ...
            geometry_params, domain_params, hydraulic_params, t_params, lognrnd_param_definer)
        num_domains = geometry_params.num_domains;
        
        mu = zeros(num_domains, 1);
        sigma = zeros(num_domains, 1);
        for idx = 1:num_domains
            [mu(idx), sigma(idx)] = lognrnd_param_definer(idx).get_params(...
                hydraulic_params.k_sat(idx) ./ geometry_params.dz, properties_array.effective_saturation(idx));
        end
        
        log_normal_params = struct('mu', [], 'sigma', []);
        log_normal_params.mu = mu;
        log_normal_params.sigma = sigma;
        
        % Calculate breakthrough
        breakthrough = hydro_1d_2domain_exchange(scale, t_idx, geometry_params, domain_params, t_params, log_normal_params);
        leachate_intercell_array = leachate_intercell_array + breakthrough;
        
        %% TODO: Alter effective saturation inside hydro
        %%
%         intercell_flux = vertcat(scale(1, :), squeeze(leachate_intercell_array(t, 1:end-1, :))) ...
%             - squeeze(leachate_intercell_array(t, :, :));
%         properties_array.effective_saturation = alter_effective_saturation(properties_array.effective_saturation, ...
%             intercell_flux, geometry_params, hydraulic_params);
    end
        
%     function new_se = alter_effective_saturation(current_se, flx, geometry_params, hydraulic_params)
%         vert_size = size(geometry_params.is_landfill_array, 1);
%         dx = reproduce(geometry_params.dx, vert_size);
%         dy = reproduce(geometry_params.dy, vert_size);
%         column_volume_array = geometry_params.dz .* geometry_params.is_landfill_array .* dx .* dy;
%         max_water_volume = (reproduce(hydraulic_params.theta_s, vert_size) - reproduce(hydraulic_params.theta_r, vert_size)) ...
%             .* column_volume_array;
%         new_se = zeros(size(current_se));
%         idx = (max_water_volume ~= 0);
%         new_se(idx) = current_se(idx) + flx(idx) ./ max_water_volume(idx);
%         % Make sure effective saturation doesn't exceed 1
%         if ~isempty(find(new_se > 1, 1, 'first'))
%             warning('Warning! Effective saturation exceeds 1');
%         end
%     end
% 
%     function res = reproduce(inp, times)
%         res = repmat(inp, [times, 1]);
%     end
    
end