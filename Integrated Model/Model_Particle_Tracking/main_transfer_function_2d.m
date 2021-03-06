function main_transfer_function_2d()

%   Integrated modelling framework for landfills
%
%   Landfill is considered as a 2-dimensional network of streamtubes of
%   different lentgh.
%
%   Features:
%     - log-normal transport function used to model travel times;
%     - variable moisture content of columns affects conductivity;
%

    clc;
    addpath('../Common/')
    tic;

    %% Preparing data for simulation:
    
    % Get spatial_params characteristics of a landfill:
    spatial_params = define_geometry();

    % Determine probability distribution parameters corresponding to defined inputs:
    lognrnd_param_definer = log_normal_params('../Common/opt_params_wt_channel_domain.mat');

    % Fluid hydraulic parameters
    hydraulic_params = lognrnd_param_definer.hydraulic_params;
    hydraulic_params.k_sat_ref = hydraulic_params.k_sat;                    % reference conductivity
    hydraulic_params.k_sat = generate_conductivities(spatial_params);       % relative conductivity compared to reference conductivity
    hydraulic_params.d = 1;                                                 % diffusion_coefficient
    
    % Generate biogeochemical properties:
    properties_array = generate_biogeochemical_properties_2d(spatial_params);
    
    % Generate precipitation data (specific, independent of area):
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
    
    progr = floor(time_params.num_intervals / 10);
    
    %% Main loop over time
    for t = 1:time_params.num_intervals
        in_flux = precipitation_in_time_vector(t) * spatial_params.dx * spatial_params.dy;
        if (precipitation_in_time_vector(t) > 0)
            [breakthrough_tmp, breakthrough_tmp_old, properties_array] = transport_lognormal(t, ...
                                                                       properties_array, ...
                                                                       in_flux, ...
                                                                       spatial_params, ...
                                                                       hydraulic_params, ...
                                                                       time_params, ...
                                                                       lognrnd_param_definer);
            leachate_out_array(t:end, :, :) = leachate_out_array(t:end, :, :) + breakthrough_tmp;
            leachate_out_array_old(t:end, :, :) = leachate_out_array_old(t:end, :, :) + breakthrough_tmp_old;
        end
        
        flx = in_flux - squeeze(leachate_out_array(t, :, :));
        properties_array.effective_saturation = ...
            alter_effective_saturation(properties_array.effective_saturation, flx, spatial_params, hydraulic_params);

        % Display progress
        if (mod(t, progr) == 0)
            disp (t / time_params.num_intervals * 100);
        end
    end

    leachate_out = squeeze(sum(sum(leachate_out_array, 3), 2));
    leachate_out_old = squeeze(sum(sum(leachate_out_array_old, 3), 2));
    
    plot(time_params.days_elapsed, leachate_out);
    hold on;
    plot(time_params.days_elapsed, leachate_out_old, 'g');
%     plot(time_params.days_elapsed, spatial_params.num_columns * precipitation_in_time_vector, 'r');
    hold off;
    
    toc;
    return
    
    function [breakthrough, breakthrough_old, properties_array] = transport_lognormal(t, properties_array, scale, ...
            spatial_params, hydraulic_params, time_params, lognrnd_param_definer)
        num_intervals = time_params.num_intervals;
        time_discretization = time_params.time_discretization;
        
        t_vector = 0 : time_discretization : time_discretization * (num_intervals - t + 1);
        
        mu = zeros(size(spatial_params.column_height_array));
        sigma = zeros(size(spatial_params.column_height_array));
        idx_calc = (spatial_params.column_height_array ~= 0);
        
        [mu(idx_calc), sigma(idx_calc)] = lognrnd_param_definer.get_params(...
            hydraulic_params.k_sat(idx_calc) ./ spatial_params.column_height_array(idx_calc), ...
            properties_array.effective_saturation(idx_calc));
        
        breakthrough_cum = zeros(num_intervals - t + 2, size(mu, 1), size(mu, 2));
        breakthrough_cum(:, (idx_calc)) = scale * log_normal_cdf(t_vector, mu(idx_calc), sigma(idx_calc))';
        breakthrough = breakthrough_cum(2:end, :, :) - breakthrough_cum(1:end-1, :, :);
        
        %% For comparison
        mu_old = zeros(size(spatial_params.column_height_array));
        sigma_old = zeros(size(spatial_params.column_height_array));
        [mu_old(idx_calc), sigma_old(idx_calc)] = lognrnd_param_definer.get_params(...
            hydraulic_params.k_sat(idx_calc) ./ spatial_params.column_height_array(idx_calc), 0.5);

        breakthrough_cum_old = zeros(num_intervals - t + 2, size(mu, 1), size(mu, 2));
        breakthrough_cum_old(:, (idx_calc)) = scale * log_normal_cdf(t_vector, mu_old(idx_calc), sigma_old(idx_calc))';
        breakthrough_old = breakthrough_cum_old(2:end, :, :) - breakthrough_cum_old(1:end-1, :, :);
    end

    function new_se = alter_effective_saturation(current_se, flx, spatial_params, hydraulic_params)
        column_volume_array = spatial_params.column_height_array * spatial_params.dx * spatial_params.dy;
        max_water_volume = (hydraulic_params.theta_s - hydraulic_params.theta_r) * column_volume_array;
        new_se = zeros(size(current_se));
        idx = (max_water_volume ~= 0);
        new_se(idx) = current_se(idx) + flx(idx) ./ max_water_volume(idx) .* spatial_params.column_height_array(idx);
    end

    function k_sat_gen = generate_conductivities(spatial_params)
        k_sat_gen = zeros(size(spatial_params.column_height_array));    % generated hydraulic conductivities
        idx = spatial_params.column_height_array > 0;                   % cells for which conductivities are to be generated
        % SET PROBABILITY DISTRIBUTION OF CONDUCTIVITIES HERE:
        randn('seed', 1);
        pdf_mean = 1e-0;
%         k_sat_gen(idx) = lognrnd(log(pdf_mean) - 1 / 8, 1, nnz(idx), 1);        % log-normal distribution
        k_sat_gen(idx) = exprnd(1e-3 / pdf_mean, nnz(idx), 1);                  % exponential distribution
%         k_sat_gen(idx) = 2 * pdf_mean * ones(nnz(idx), 1);                      % uniform distribution
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