function main_solute_particle_tracking()
    clc;
    addpath('../Common/');
    
    tic;
    
    rng(10);
%     rand('seed', 10);
    
    %% Define inputs
    % Solute transport params
    sol = struct();
    sol.u = 5e-1;           % flow rate
    sol.d = 1e-1;           % diffusion parameter
    sol.dv = 1e-4;          % solute unit volume (particle volume)

    % Time params
    dt = 0.01;
    t_end = 10;
    
    % Geometry
    geom = struct();
    geom.zn = 10;
    geom.dx = 1;                                % cell length
    geom.dy = 1;                                % cell width
    geom.dz = 1;                                % cell height
    geom.dv = geom.dx * geom.dy * geom.dz;      % cell volume
    geom.l = geom.dz * ones(geom.zn, 1);         % vector of all cells' heights
    geom.c_ini = 0.03 * ones(geom.zn, 1);         % initial concentrations of solute
    geom.is_landfill_array = ones(geom.zn, 1);   % vector denoting whether the cell belongs to landfill, or not
    
    % Precipitation and time parameters
    file_name = '../Data/precipitation_daily_KNMI_20110908.txt';
    [precipitation_intensity_time_vector, time_params, start_date] = read_precipitation_data_csv(file_name);
    precipitation_intensity_time_vector = zeros(size(precipitation_intensity_time_vector));
    precipitation_intensity_time_vector(1) = 1e-3;
%     t = time_params.t;
%     dt = time_params.time_discretization;

    % Determine probability distribution parameters corresponding to defined inputs:
    lognrnd_param_definer = log_normal_params('../Common/opt_params_wt_channel_domain.mat');
    
    % Fluid hydraulic parameters
    hydraulic_params = lognrnd_param_definer.hydraulic_params;
    hydraulic_params.k_sat_ref = hydraulic_params.k_sat;    % reference conductivity
    hydraulic_params.k_sat = 5e-2;                            % relative conductivity compared to reference conductivity
    hydraulic_params.d = 1;                                 % diffusion_coefficient

    % initial SE
    properties_array = generate_biogeochemical_properties_3d(geom, hydraulic_params);
    
    %% Initialize result
    leachate_flux = zeros(time_params.num_intervals, numel(geom.is_landfill_array));
    solute_out_flux = zeros(time_params.num_intervals, 1);

    %% Precompute 
    particles_arr = solute_particle_class(geom, sol);

    % Progress step
    progr = floor(time_params.num_intervals / 10);
    
    %% Main loop
    for t_idx = 1:time_params.num_intervals
        [leachate_flux, properties_array] = transport_lognormal(leachate_flux, t_idx, properties_array, ...
            precipitation_intensity_time_vector(t_idx), geom, hydraulic_params, time_params, lognrnd_param_definer);
        
        % Display progress
        if (mod(t_idx, progr) == 0)
            disp (t_idx / time_params.num_intervals * 100);
        end
    end

    solute_out_flux = solute_out_flux * 1e-1;
    
    %% Display results
    plot(time_params.days_elapsed, leachate_flux(:, end), 'b');
    hold on;
    sof = smooth(solute_out_flux * sol.dv, 1);
    plot(time_params.days_elapsed, 1e-0 * sof, 'r');
    plot(time_params.days_elapsed, 1e-0 * sof .* sol.dv ./ (leachate_flux(:, end) + sof .* sol.dv), 'k');
    hold off;
    
    toc;
    
    return
    
    function [leachate_intercell_array, properties_array] = transport_lognormal(leachate_intercell_array, t_idx, properties_array, scale, ...
            spatial_params, hydraulic_params, time_params, lognrnd_param_definer)
        [mu, sigma] = lognrnd_param_definer.get_params(...
            hydraulic_params.k_sat ./ spatial_params.dz, properties_array.effective_saturation);
        skip_cell_idx = (spatial_params.is_landfill_array == 0);
        mu(skip_cell_idx) = 0;
        sigma(skip_cell_idx) = Inf;
        
        t_vector = 0 : time_params.time_discretization : time_params.time_discretization * (time_params.num_intervals - t_idx + 1);
        
        scale = scale * spatial_params.dx * spatial_params.dy;
        breakthrough_cum = scale * log_normal_cdf(t_vector, mu(1), sigma(1));
        breakthrough = breakthrough_cum(2:end) - breakthrough_cum(1:end-1);
        leachate_intercell_array(t_idx:end, 1) = leachate_intercell_array(t_idx:end, 1) + breakthrough';
        for cell_idx = 2:spatial_params.zn
            breakthrough_cum = leachate_intercell_array(t_idx, cell_idx - 1) ...
                * log_normal_cdf(t_vector, mu(cell_idx), sigma(cell_idx));
            breakthrough = breakthrough_cum(2:end) - breakthrough_cum(1:end-1);
            leachate_intercell_array(t_idx:end, cell_idx) = leachate_intercell_array(t_idx:end, cell_idx) + breakthrough';
        end
        
        intercell_flux = horzcat(scale, leachate_intercell_array(t_idx, 1:end-1)) - leachate_intercell_array(t_idx, :);
        properties_array.effective_saturation = alter_effective_saturation(properties_array.effective_saturation, ...
            intercell_flux, spatial_params, hydraulic_params);
        
        %% Solute transport module:
%         % Average water velocity per cell - depends on conductivity:
%         mean_water_velocity = spatial_params.l ./ exp(mu + sigma .^ 2 / 2);
%         mean_water_velocity(mean_water_velocity == Inf) = 0;
        % Actual velocity is assumed to be ###
        water_velocity = intercell_flux' ./ (spatial_params.dx * spatial_params.dy) ./ properties_array.effective_saturation;
        
        [particles_arr, solute_out_flux(t_idx)] = particles_arr.update(time_params.t(t_idx), water_velocity);
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