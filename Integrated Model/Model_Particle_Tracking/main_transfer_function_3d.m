function main_transfer_function_3d()

%   Integrated modelling framework for landfills
%
%   Landfill is considered as a 3-dimensional network of cells of the same
%   volume. Water flows through these cells.
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
    properties_array = generate_biogeochemical_properties_3d(spatial_params);
    
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
	leachate_intercell_array = zeros(time_params.num_intervals, spatial_params.xn, spatial_params.yn, spatial_params.zn);
    leachate_out = zeros(1, time_params.num_intervals);
    
    %% Main loop over time
    for t = 1:time_params.num_intervals
        in_flux = precipitation_in_time_vector(t) .* (spatial_params.column_height_array > 0); % input flux per column
        [leachate_intercell_array(t:end, :, :, :), properties_array] = ...
            transport_lognormal(leachate_intercell_array(t:end, :, :, :), ...
            t, ...
            properties_array, ...
            in_flux, ...
            spatial_params, ...
            hydraulic_params, ...
            time_params, ...
            lognrnd_param_definer);
        
        tmp_dims = size(leachate_intercell_array(t, :, :, :));
        intercell_flux = cat(3, in_flux, reshape(leachate_intercell_array(t, :, :, 1:end-1), [tmp_dims(2:3) (tmp_dims(4)-1)])) - ...
            reshape(leachate_intercell_array(t, :, :, :), tmp_dims(2:4));
        properties_array.effective_saturation = ...
            alter_effective_saturation(properties_array.effective_saturation, intercell_flux, spatial_params, hydraulic_params);

        % Display progress
        if (mod(t, 1000) == 0)
            disp (t / time_params.num_intervals * 100);
        end
    end

    leachate_out = squeeze(sum(sum(leachate_intercell_array(:, :, :, end), 3), 2));
    hold on;
    plot(time_params.days_elapsed, leachate_out, 'r');
    hold off;
    
    return
    
    function [leachate_intercell_array, properties_array] = transport_lognormal(leachate_intercell_array, t, properties_array, scale, ...
            spatial_params, hydraulic_params, time_params, lognrnd_param_definer)
        num_intervals = time_params.num_intervals;
        time_discretization = time_params.time_discretization;
        t_vector = 0 : time_discretization : time_discretization * (num_intervals - t);
        scale = repmat(scale, [1, 1, numel(t_vector)]);
        scale = permute(scale, [3, 1, 2]);
        
        mu = zeros(size(spatial_params.is_landfill_array));
        sigma = zeros(size(spatial_params.is_landfill_array));
        idx_calc = (spatial_params.column_height_array ~= 0);
        idx_calc_3d = repmat(idx_calc, [1, 1, spatial_params.zn]);
        
        [mu(idx_calc_3d), sigma(idx_calc_3d)] = lognrnd_param_definer.get_params(...
            hydraulic_params.k_sat ./ spatial_params.dz, ...
            properties_array.effective_saturation(idx_calc_3d));
        
        % Initialize breakthrough but leave breakthrough(1, :, :) = 0;
        breakthrough = zeros(num_intervals - t + 1, size(mu, 1), size(mu, 2));
        
        t_idx = 2:num_intervals - t + 1;
        % Calculations for the upper layer:
        idx_calc_layer = logical(zeros(size(mu)));
        idx_calc_layer(:, :, 1) = idx_calc;
        breakthrough(t_idx, idx_calc) = scale(t_idx, idx_calc) .* time_discretization .* ...
            log_normal_pdf(t_vector(t_idx), mu(idx_calc_layer), sigma(idx_calc_layer))';
        breakthrough = avg_flow(breakthrough);
        leachate_intercell_array(:, :, :, 1) = leachate_intercell_array(:, :, :, 1) + breakthrough;
        % Calculations for other layers:
        for layer_idx = 2:spatial_params.zn
            idx_calc_layer = logical(zeros(size(mu)));
            idx_calc_layer(:, :, layer_idx) = idx_calc;
            scale = leachate_intercell_array(1, :, :, layer_idx-1);
            scale = repmat(scale, [1, 1, numel(t_vector)]);
            scale = permute(scale, [3, 1, 2]);
            breakthrough(t_idx, idx_calc) = scale(t_idx, idx_calc) .* time_discretization .* ...
                log_normal_pdf(t_vector(t_idx), mu(idx_calc_layer), sigma(idx_calc_layer))';
            breakthrough = avg_flow(breakthrough);
            leachate_intercell_array(:, :, :, layer_idx) = leachate_intercell_array(:, :, :, layer_idx) + breakthrough;
        end
    end

    function new_se = alter_effective_saturation(current_se, flx, spatial_params, hydraulic_params)
        cell_volume = spatial_params.dx * spatial_params.dy * spatial_params.dz .* spatial_params.is_landfill_array;
        max_water_volume = (hydraulic_params.theta_s - hydraulic_params.theta_r) * cell_volume;
        new_se = zeros(size(current_se));
        idx = (max_water_volume ~= 0);
        new_se(idx) = current_se(idx) + flx(idx) ./ max_water_volume(idx);
    end

%     %% TODO:
%     function res = solute_transport(properties_array, spatial_params, hydraulic_params, time_params)
%         
%         d = hydraulic_params.d;
%         dt = time_params.time_discretization;
%         l = spatial_params.column_height_array;
%         v = properties_array.mean_v;
%         
%         
%         
%         solute_concentration = -(1 ./ 2) .* erf((1 ./ 2) .* (-x + v .* dt) ./ (sqrt(dt) .* sqrt(d))) + ...
%             (1 ./ 2) .* erf((1 ./ 2) .* (-x + l + v .* dt) ./ (sqrt(dt) .* sqrt(d)));
%     end

end