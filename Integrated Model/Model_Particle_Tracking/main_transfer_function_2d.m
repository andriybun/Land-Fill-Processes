function main_transfer_function_2d()

%   Integrated modelling framework for landfills
%
%   Features:
%     - log-normal transport function used to model travel times;
%     - variable moisture content of columns affects conductivity;

%% TODO: warning if time step is too big and there is a danger of overflow

    clc;
    addpath('../Common/')

    %% Input parameters:
    
    % Time
    start_date = struct();
    start_date.year = 2012;
    start_date.month = 1;
    start_date.day = 1;
    time_params.max_days = 30; % 365;                                                           % number of simulation days
    time_params.time_discretization = 3600;                                                     % in seconds
    time_params.intervals_per_day = 24 * 3600 / time_params.time_discretization;
    num_intervals = time_params.max_days * time_params.intervals_per_day;                       % in {time step}
    time_params.num_intervals = num_intervals;
    time_params.days_elapsed = (0 : 1: (num_intervals-1)) / time_params.intervals_per_day;
    
    % Fluid velocity parameters
    hydraulic_params.k_sat = 1e-4;                      % m / s
    hydraulic_params.theta_r = 0.102;                   % residual water content
    hydraulic_params.theta_s = 0.368;                   % saturated water content
    hydraulic_params.d = 1;                             % diffusion_coefficient
%     expected_fluid_velocity = ...
%         expected_fluid_velocity_mps * time_discretization;  % m / {time step}
%     variance = 3 * expected_fluid_velocity;                 % 95% confidence interval width

    %% Preparing data for simulation:
    
    % Get spatial_params characteristics of a landfill:
    spatial_params = define_geometry();

    % Determine probability distribution parameters corresponding to defined inputs:
    lognrnd_param_definer = log_normal_params();

    % Generate biogeochemical properties:
    properties_array = generate_biogeochemical_properties_2d(spatial_params);
    
    % Generate precipitation data (specific, independent of area):
    rand('seed', 1);
%     precipitation_intensity_time_vector = generate_precipitation_data(start_date, time_params);
    file_name = '../Data/precipitation_daily_KNMI_20110908.txt';
    [precipitation_intensity_time_vector, time_params, start_date] = read_precipitation_data_csv(file_name);
%     file_name = '../Data/precip_braambergen.mat';
%     [precipitation_intensity_time_vector, time_params, start_date] = read_precipitation_data_braambergen(file_name);
    
    % Net amounts of precipitation from rainfall events entering each column / pathway:
    precipitation_in_time_vector = precipitation_intensity_time_vector * ...
        (spatial_params.dx * spatial_params.dy);
    
    % Array specifying amounts of leachate leaving the landfill
	leachate_out_array = zeros(time_params.num_intervals, spatial_params.xn, spatial_params.yn);
    leachate_out_array_old = zeros(time_params.num_intervals, spatial_params.xn, spatial_params.yn);
    leachate_out = zeros(1, time_params.num_intervals);
    
    %% Main loop over time
    for t = 1:time_params.num_intervals
        in_flux = precipitation_in_time_vector(t); % * time_params.time_discretization;
        if (precipitation_in_time_vector(t) > 0)
            [breakghrough_tmp, breakthrough_tmp_old, properties_array] = transport_lognormal(t, ...
                                                                       properties_array, ...
                                                                       in_flux, ...
                                                                       spatial_params, ...
                                                                       hydraulic_params, ...
                                                                       time_params, ...
                                                                       lognrnd_param_definer);
            leachate_out_array(t:end, :, :) = leachate_out_array(t:end, :, :) + breakghrough_tmp;
            leachate_out_array_old(t:end, :, :) = leachate_out_array_old(t:end, :, :) + breakthrough_tmp_old;
        end
        
        flx = in_flux - squeeze(leachate_out_array(t, :, :));
        properties_array.effective_saturation = ...
            alter_effective_saturation(properties_array.effective_saturation, flx, spatial_params, hydraulic_params);

        % Display progress
        if (mod(t, 1000) == 0)
            disp (t / time_params.num_intervals * 100);
        end
    end

    leachate_out = squeeze(sum(sum(leachate_out_array, 3), 2));
    leachate_out_old = squeeze(sum(sum(leachate_out_array_old, 3), 2));
    
    hold on;
    plot(time_params.days_elapsed, leachate_out);
    plot(time_params.days_elapsed, leachate_out_old, 'g');
%     plot(time_params.days_elapsed, spatial_params.num_columns * precipitation_in_time_vector, 'r');
    hold off;
    
    return
    
    function [breakthrough, breakthrough_old, properties_array] = transport_lognormal(t, properties_array, scale, ...
            spatial_params, hydraulic_params, time_params, lognrnd_param_definer)
        num_intervals = time_params.num_intervals;
        time_discretization = time_params.time_discretization;
        
        t_vector = 0 : time_discretization : time_discretization * (num_intervals - t);
        
        mu = zeros(size(spatial_params.column_height_array));
        sigma = zeros(size(spatial_params.column_height_array));
        idx_calc = (spatial_params.column_height_array ~= 0);
        
        [mu(idx_calc), sigma(idx_calc)] = lognrnd_param_definer.get_params(...
            hydraulic_params.k_sat ./ spatial_params.column_height_array(idx_calc), ...
            properties_array.effective_saturation(idx_calc));
        breakthrough = zeros(num_intervals - t + 1, size(mu, 1), size(mu, 2));
        % Leave breakthrough(1, :, :) = 0;
        t_idx = 2:num_intervals - t + 1;
        breakthrough(t_idx, (idx_calc)) = scale * time_discretization * log_normal_pdf(t_vector(t_idx), mu(idx_calc), sigma(idx_calc))';
%         for t_idx = 2:num_intervals - t + 1
%             breakthrough(t_idx, (idx_calc)) = scale * time_discretization * exp(-(log(t_vector(t_idx)) - mu(idx_calc)).^2 ./ ...
%                 (2 .* sigma(idx_calc) .* sigma(idx_calc))) ./ (sqrt(2 * pi) .* sigma(idx_calc) .* t_vector(t_idx));
%         end
        
        %% For comparison
        mu_old = zeros(size(spatial_params.column_height_array));
        sigma_old = zeros(size(spatial_params.column_height_array));
        [mu_old(idx_calc), sigma_old(idx_calc)] = lognrnd_param_definer.get_params(...
            hydraulic_params.k_sat ./ spatial_params.column_height_array(idx_calc), 0.7);
        breakthrough_old = zeros(num_intervals - t + 1, size(mu_old, 1), size(mu_old, 2));
        % Leave breakthrough_old(1, :, :) = 0;
        breakthrough_old(t_idx, (idx_calc)) = scale * time_discretization * log_normal_pdf(t_vector(t_idx), mu_old(idx_calc), sigma_old(idx_calc))';
%         for t_idx = 2:num_intervals - t + 1
%             breakthrough_old(t_idx, (idx_calc)) = scale * time_discretization * exp(-(log(t_vector(t_idx)) - mu_old(idx_calc)).^2 ./ ...
%                 (2 .* sigma_old(idx_calc) .* sigma_old(idx_calc))) ./ (sqrt(2 * pi) .* sigma_old(idx_calc) .* t_vector(t_idx));
%         end
        
    end

    function new_se = alter_effective_saturation(current_se, flx, spatial_params, hydraulic_params)
        column_volume_array = spatial_params.column_height_array * spatial_params.dx * spatial_params.dy;
        max_water_volume = (hydraulic_params.theta_s - hydraulic_params.theta_r) * column_volume_array;
        new_se = zeros(size(current_se));
        idx = (max_water_volume ~= 0);
        new_se(idx) = current_se(idx) + flx(idx) ./ max_water_volume(idx);
    end

    %% TODO:
    function res = solute_transport(properties_array, spatial_params, hydraulic_params, time_params)
        
        d = hydraulic_params.d;
        dt = time_params.time_discretization;
        l = spatial_params.column_height_array;
        v = properties_array.mean_v;
        
        
        
        solute_concentration = -(1 ./ 2) .* erf((1 ./ 2) .* (-x + v .* dt) ./ (sqrt(dt) .* sqrt(d))) + ...
            (1 ./ 2) .* erf((1 ./ 2) .* (-x + l + v .* dt) ./ (sqrt(dt) .* sqrt(d)));
    end

end