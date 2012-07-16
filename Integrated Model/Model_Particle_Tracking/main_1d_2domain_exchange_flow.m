function main_1d_2domain_exchange_flow()
    clc;
    tic;
    addpath('../Common/');

    % Time params
    start_date = struct();
    start_date.year = 2012;
    start_date.month = 1;
    start_date.day = 1;
    time_params.max_days = 100;                                                           % number of simulation days
    time_params.intervals_per_day = 24;
    time_params.time_discretization = 1 / time_params.intervals_per_day;                 % in days
    time_params.num_intervals = time_params.max_days * time_params.intervals_per_day;    % in {time step}
    time_params.days_elapsed = (0 : 1: (time_params.num_intervals-1)) / time_params.intervals_per_day;
    num_intervals = time_params.num_intervals;
    
    % Geometry params
    spatial_params = struct();
    spatial_params.dx = [1, 1];
    chan_width = 0.1;
    spatial_params.dy = [1 - chan_width, chan_width];
    spatial_params.dz = 1;
    spatial_params.zn = 10;
    spatial_params.is_landfill_array = ones(spatial_params.zn, 2);
    
    % Characteristics of eschange between domains (rate of flow from matrix
    % domain to channel domain)
    z = cumsum(ones(size(spatial_params.is_landfill_array(:, 1))), 1) - spatial_params.dz;
    exchange_rate = 0.3 * (1 - 1 / max(z)^(1/2) * z.^(1/2));
%     exchange_rate = 0.3 * (1 - 1 / max(z) * z);
    
    file_name = '../Data/precipitation_daily_KNMI_20110908.txt';
    [precipitation_intensity_time_vector, time_params, start_date] = read_precipitation_data_csv(file_name);
    num_intervals = time_params.num_intervals;
    %% TESTING (let all the water flow out to check mass balance):
    extra_days = 100;
    time_params.max_days = time_params.max_days + extra_days;
    time_params.num_intervals = time_params.max_days * time_params.intervals_per_day;                           % in {time step}
    num_intervals = time_params.num_intervals;
    time_params.days_elapsed = (0 : 1: (time_params.num_intervals-1)) / time_params.intervals_per_day;
    precipitation_intensity_time_vector = cat(2, precipitation_intensity_time_vector, zeros(1, extra_days * time_params.intervals_per_day));
    %% END TESTING
    
    precipitation_intensity_time_vector = zeros(1, num_intervals);
    precipitation_intensity_time_vector(1) = 1e-4;
    
    % Determine probability distribution parameters corresponding to defined inputs:
    lognrnd_param_definer(1) = log_normal_params('../Common/opt_params_wt_matrix_domain.mat');
    lognrnd_param_definer(2) = log_normal_params('../Common/opt_params_wt_channel_domain.mat');
    
    % Fluid velocity parameters (1st - matrix domain; 2nd - channel domain)
    % Note: vectors must be row vectors
    hydraulic_params            = aos_to_soa([lognrnd_param_definer(1).hydraulic_params
                                              lognrnd_param_definer(2).hydraulic_params]');
    hydraulic_params.k_sat_ref  = hydraulic_params.k_sat;
    hydraulic_params.k_sat      = [50, 1e-2];
	hydraulic_params.d          = [1, 1];                                                        % diffusion_coefficient

    % initial SE
    properties_array = generate_biogeochemical_properties_3d(spatial_params, hydraulic_params);
    
    % Result
    leachate_flux = zeros(horzcat(num_intervals, size(spatial_params.is_landfill_array)));
    
    progr = floor(time_params.num_intervals / 10);
    
    for t = 1:num_intervals
        [leachate_flux, properties_array] = transport_lognormal(leachate_flux, t, ...
            properties_array, precipitation_intensity_time_vector(t), spatial_params, ...
            hydraulic_params, time_params, lognrnd_param_definer);
        
        % Display progress
        if (mod(t, progr) == 0)
            disp (t / time_params.num_intervals * 100);
        end
    end

%     plot(time_params.days_elapsed, sum(leachate_flux, 3));
    plot(time_params.days_elapsed, squeeze(leachate_flux(:, end, 1)), 'r');
    hold on;
    plot(time_params.days_elapsed, squeeze(leachate_flux(:, end, 2)), 'g');
    plot(time_params.days_elapsed, squeeze(sum(leachate_flux(:, end, :), 3)), 'k');
    hold off;
%     hold on;
%     plot(time_params.days_elapsed, leachate_flux_matrix_domain, 'g');
%     plot(time_params.days_elapsed, leachate_flux_fast_domain, 'r');
%     hold off;
%     figure(2);
%     mesh(se);

    toc;

    return 
    
    function [leachate_intercell_array, properties_array] = transport_lognormal(leachate_intercell_array, t, properties_array, scale, ...
            spatial_params, hydraulic_params, time_params, lognrnd_param_definer)
        % Number of domains
        num_domains = numel(lognrnd_param_definer);
        
        % Define log-normal parameters
        mu = zeros(size(properties_array.effective_saturation));
        sigma = zeros(size(properties_array.effective_saturation));
        for idx = 1:num_domains
            [mu(:, idx), sigma(:, idx)] = lognrnd_param_definer(idx).get_params(...
                hydraulic_params.k_sat(idx) ./ spatial_params.dz, properties_array.effective_saturation(:, idx));
        end
        skip_cell_idx = (spatial_params.dz == 0);
        mu(skip_cell_idx) = 0;
        sigma(skip_cell_idx) = Inf;
        
        % Time vector
        t_vector = 0 : time_params.time_discretization : time_params.time_discretization * (num_intervals - t + 1);
        
        % Water input at time t
        scale = scale .* spatial_params.dx .* spatial_params.dy;
        exchange = exchange_rate(1) * scale(1);
        scale = scale + [-exchange, exchange];
        scale = reproduce(scale, numel(t_vector));
        
        % Calculate breakthrough (upper layer)
        breakthrough_cum = scale .* log_normal_cdf(t_vector, mu(1, :)', sigma(1, :)')';
        breakthrough = breakthrough_cum(2:end, :) - breakthrough_cum(1:end-1, :);
        leachate_intercell_array(t:end, 1, :) = leachate_intercell_array(t:end, 1, :) ...
            + permute(breakthrough, [1 3 2]);
        
        % Calculate breakthrough (other layers)
        for cellIdx = 2:spatial_params.zn
            flux = squeeze(leachate_intercell_array(t, cellIdx - 1, :))';
            exchange = exchange_rate(cellIdx) * flux(1);
            flux = flux + [-exchange, exchange];
            flux = reproduce(flux, numel(t_vector));
            breakthrough_cum = flux .* log_normal_cdf(t_vector, mu(cellIdx, :)', sigma(cellIdx, :)')';
            breakthrough = breakthrough_cum(2:end, :) - breakthrough_cum(1:end-1, :);
            leachate_intercell_array(t:end, cellIdx, :) = leachate_intercell_array(t:end, cellIdx, :) ...
                + permute(breakthrough, [1 3 2]);
        end
        
        intercell_flux = vertcat(scale(1, :), squeeze(leachate_intercell_array(t, 1:end-1, :))) ...
            - squeeze(leachate_intercell_array(t, :, :));
        properties_array.effective_saturation = alter_effective_saturation(properties_array.effective_saturation, ...
            intercell_flux, spatial_params, hydraulic_params);
    end
        
    function new_se = alter_effective_saturation(current_se, flx, spatial_params, hydraulic_params)
        vert_size = size(spatial_params.is_landfill_array, 1);
        dx = reproduce(spatial_params.dx, vert_size);
        dy = reproduce(spatial_params.dy, vert_size);
        column_volume_array = spatial_params.dz .* spatial_params.is_landfill_array .* dx .* dy;
        max_water_volume = (reproduce(hydraulic_params.theta_s, vert_size) - reproduce(hydraulic_params.theta_r, vert_size)) ...
            .* column_volume_array;
        new_se = zeros(size(current_se));
        idx = (max_water_volume ~= 0);
        new_se(idx) = current_se(idx) + flx(idx) ./ max_water_volume(idx);
        % Make sure effective saturation doesn't exceed 1
        if ~isempty(find(new_se > 1, 1, 'first'))
            warning('Warning! Effective saturation exceeds 1');
        end
    end

    function res = reproduce(inp, times)
        res = repmat(inp, [times, 1]);
    end
    
end