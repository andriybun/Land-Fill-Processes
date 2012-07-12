function try_1dim_flow()
    clc;
    addpath('../Common/');

    % Time params
    start_date = struct();
    start_date.year = 2012;
    start_date.month = 1;
    start_date.day = 1;
    time_params.max_days = 30;                                                           % number of simulation days
    time_params.time_discretization = 3600;                                                     % in seconds
    time_params.intervals_per_day = 24 * 3600 / time_params.time_discretization;
    num_intervals = time_params.max_days * time_params.intervals_per_day;                       % in {time step}
    time_params.num_intervals = num_intervals;
    time_params.days_elapsed = (0 : 1: (num_intervals-1)) / time_params.intervals_per_day;
    
    % Flow params
    spatial_params = struct();
    spatial_params.dx = 1;
    spatial_params.dy = 1;
    spatial_params.dz = ones(8, 1);                     % vertical dimensions of spatial_params.dz, forming column
    spatial_params.dz(5:end) = 0;
    
    properties_array = struct();
    properties_array.effective_saturation = 0.3 * ones(size(spatial_params.dz));    % initial SE
    
%     file_name = '../Data/precipitation_daily_KNMI_20110908.txt';
%     [precipitation_intensity_time_vector, time_params, start_date] = read_precipitation_data_csv(file_name);
    precipitation_intensity_time_vector = zeros(1, num_intervals);
    precipitation_intensity_time_vector(1) = 1e-4;
    
    % Determine probability distribution parameters corresponding to defined inputs:
    lognrnd_param_definer = log_normal_params('../Common/opt_params_wt_matrix_domain.mat');
    
    % Fluid hydraulic parameters
    hydraulic_params = lognrnd_param_definer.hydraulic_params;
    hydraulic_params.k_sat_ref = hydraulic_params.k_sat;    % reference conductivity
    hydraulic_params.k_sat = 1;                             % relative conductivity compared to reference conductivity
    hydraulic_params.d = 1;                                 % diffusion_coefficient
    
    % Result
    leachate_flux = zeros(num_intervals, numel(spatial_params.dz));
%     se = zeros(numel(spatial_params.dz), num_intervals);
    
    for t = 1:num_intervals
        [leachate_flux, properties_array] = transport_lognormal(leachate_flux, t, properties_array, ...
            precipitation_intensity_time_vector(t), spatial_params, hydraulic_params, time_params, lognrnd_param_definer);
%         se(:, t) = effective_saturation;
    end

    plot(time_params.days_elapsed, leachate_flux);
%     figure(2);
%     mesh(se);

    return 
    
    function [leachate_intercell_array, properties_array] = transport_lognormal(leachate_intercell_array, t, properties_array, scale, ...
            spatial_params, hydraulic_params, time_params, lognrnd_param_definer)
        [mu, sigma] = lognrnd_param_definer.get_params(...
            hydraulic_params.k_sat ./ spatial_params.dz, properties_array.effective_saturation);
        skip_cell_idx = (spatial_params.dz == 0);
        mu(skip_cell_idx) = 0;
        sigma(skip_cell_idx) = Inf;
        
        t_vector = 0 : time_params.time_discretization : time_params.time_discretization * (num_intervals - t + 1);
        
        scale = scale * spatial_params.dx * spatial_params.dy;
        breakthrough_cum = scale * log_normal_cdf(t_vector, mu(1), sigma(1));
        breakthrough = breakthrough_cum(2:end) - breakthrough_cum(1:end-1);
        leachate_intercell_array(t:end, 1) = leachate_intercell_array(t:end, 1) + breakthrough';
        for cell_idx = 2:numel(spatial_params.dz)
            breakthrough_cum = leachate_intercell_array(t, cell_idx - 1) ...
                * log_normal_cdf(t_vector, mu(cell_idx), sigma(cell_idx));
            breakthrough = breakthrough_cum(2:end) - breakthrough_cum(1:end-1);
            leachate_intercell_array(t:end, cell_idx) = leachate_intercell_array(t:end, cell_idx) + breakthrough';
        end
        
        intercell_flux = horzcat(scale, leachate_intercell_array(t, 1:end-1)) - leachate_intercell_array(t, :);
        properties_array.effective_saturation = alter_effective_saturation(properties_array.effective_saturation, ...
            intercell_flux, spatial_params, hydraulic_params);
    end
        
    function new_se = alter_effective_saturation(current_se, flx, spatial_params, hydraulic_params)
        column_volume_array = spatial_params.dz .* spatial_params.dx .* spatial_params.dy;
        max_water_volume = (hydraulic_params.theta_s - hydraulic_params.theta_r) * column_volume_array;
        new_se = zeros(size(current_se));
        idx = (max_water_volume ~= 0);
        new_se(idx) = current_se(idx) + flx(idx)' ./ max_water_volume(idx);
        % Make sure effective saturation doesn't exceed 1
        if ~isempty(find(new_se > 1, 1, 'first'))
            warning('Warning! Effective saturation exceeds 1');
        end
    end
    
end