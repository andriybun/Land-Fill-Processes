function main_1d_flow()
    clc;
    addpath('../Common/');

    % Time params
    start_date = struct();
    start_date.year = 2012;
    start_date.month = 1;
    start_date.day = 1;
    time_params.max_days = 250;                                                         % number of simulation days
    time_params.intervals_per_day = 24;
    time_params.time_discretization = 1 / time_params.intervals_per_day;                 % in days
    time_params.num_intervals = time_params.max_days * time_params.intervals_per_day;    % in {time step}
    time_params.days_elapsed = (0 : 1: (time_params.num_intervals-1)) / time_params.intervals_per_day;
    num_intervals = time_params.num_intervals;
    
    % Flow params
    spatial_params = struct();
    spatial_params.dx = 1;
    spatial_params.dy = 1;
    spatial_params.dz = 1;                              % vertical dimensions of spatial_params.dz, forming column
    spatial_params.zn = 8;
    spatial_params.is_landfill_array = ones(spatial_params.zn, 1);
    spatial_params.is_landfill_array(5:end) = 0;
    
%     file_name = '../Data/precipitation_daily_KNMI_20110908.txt';
%     [precipitation_intensity_time_vector, time_params, start_date] = read_precipitation_data_csv(file_name);
    precipitation_intensity_time_vector = zeros(1, num_intervals);
    precipitation_intensity_time_vector(1) = 1e-4;
    
    % Determine probability distribution parameters corresponding to defined inputs:
    lognrnd_param_definer = log_normal_params('../Common/opt_params_wt_matrix_domain.mat');
    
    % Fluid hydraulic parameters
    hydraulic_params = lognrnd_param_definer.hydraulic_params;
    hydraulic_params.k_sat_ref = hydraulic_params.k_sat;    % reference conductivity
    hydraulic_params.k_sat = 50;                            % relative conductivity compared to reference conductivity
    hydraulic_params.d = 1;                                 % diffusion_coefficient

	% initial SE
    properties_array = generate_biogeochemical_properties_3d(spatial_params, hydraulic_params);

    % Result
    leachate_flux = zeros(num_intervals, numel(spatial_params.is_landfill_array));

    progr = floor(time_params.num_intervals / 10);
    
    for t = 1:num_intervals
        [leachate_flux, properties_array] = transport_lognormal(leachate_flux, t, properties_array, ...
            precipitation_intensity_time_vector(t), spatial_params, hydraulic_params, time_params, lognrnd_param_definer);
        % Display progress
        if (mod(t, progr) == 0)
            disp (t / time_params.num_intervals * 100);
        end
    end

    plot(time_params.days_elapsed, leachate_flux);

    return 
    
    function [leachate_intercell_array, properties_array] = transport_lognormal(leachate_intercell_array, t, properties_array, scale, ...
            spatial_params, hydraulic_params, time_params, lognrnd_param_definer)
        [mu, sigma] = lognrnd_param_definer.get_params(...
            hydraulic_params.k_sat ./ spatial_params.dz, properties_array.effective_saturation);
        skip_cell_idx = (spatial_params.is_landfill_array == 0);
        mu(skip_cell_idx) = 0;
        sigma(skip_cell_idx) = Inf;
        
        t_vector = 0 : time_params.time_discretization : time_params.time_discretization * (num_intervals - t + 1);
        
        scale = scale * spatial_params.dx * spatial_params.dy;
        breakthrough_cum = scale * log_normal_cdf(t_vector, mu(1), sigma(1));
        breakthrough = breakthrough_cum(2:end) - breakthrough_cum(1:end-1);
        leachate_intercell_array(t:end, 1) = leachate_intercell_array(t:end, 1) + breakthrough';
        for cell_idx = 2:spatial_params.zn
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
        column_volume_array = spatial_params.dz .* spatial_params.is_landfill_array .* spatial_params.dx .* spatial_params.dy;
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