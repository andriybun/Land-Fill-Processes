function main_1d_flow()
    clc;
    addpath('../Common/');
    rng(1);

    tic;
    
    % Flow params
    geom_params = struct();
    geom_params.dx = 1;
    geom_params.dy = 1;
    geom_params.dz = 1;                              % vertical dimensions of geom_params.dz, forming column
    geom_params.xn = 1;
    geom_params.yn = 1;
    geom_params.zn = 10;
    if geom_params.yn == 1 && geom_params.xn > 1
        error('Dimensions error: yn cannot be one while xn is greater than 1. Try switching values.');
    end
    geom_params.is_landfill_array = ones(geom_params.zn, geom_params.yn, geom_params.xn);
%     geom_params.is_landfill_array(5:end) = 0;
    


    % Precipitation data
    file_name = '../Data/precipitation_daily_KNMI_20110908.txt';
    [precipitation_intensity_time_vector, time_params, start_date] = read_precipitation_data_csv(file_name);
    t = time_params.t;
%     precipitation_intensity_time_vector = zeros(1, time_params.num_intervals);
%     precipitation_intensity_time_vector(1) = 1e-4;
    
    % Determine probability distribution parameters corresponding to defined inputs:
    lognrnd_param_definer = log_normal_params('../Common/opt_params_wt_channel_domain.mat');
    
    % Hydraulic parameters
    hydraulic_params = lognrnd_param_definer.hydraulic_params;
    hydraulic_params.k_sat_ref = hydraulic_params.k_sat;    % reference conductivity
    hydraulic_params.k_sat = 1e-2;                          % relative conductivity compared to reference conductivity
    hydraulic_params.d = 1;                                 % diffusion_coefficient    
    
    % Bio-chemical composition parameters
    properties_array = generate_biogeochemical_properties_3d(geom_params, hydraulic_params);
    effective_saturation = properties_array.effective_saturation;
% es = zeros(geom_params.zn, time_params.num_intervals);
% es(:, 1) = effective_saturation;
    %% 
    dv = 1e-5;
    wp = water_particle_class(geom_params, dv, hydraulic_params, lognrnd_param_definer);
    
    % Result
    leachate_flux = zeros(1, time_params.num_intervals);

    progr = floor(time_params.num_intervals / 10);
    for t_idx = 1:time_params.num_intervals
        % Influx per cell
        in_flux = precipitation_intensity_time_vector(t_idx) * geom_params.dy * geom_params.dx;
        [wp, leachate_flux(t_idx)] = wp.add_water_in(t(t_idx), in_flux, effective_saturation);
        conc = wp.get_concentrations();
        effective_saturation = alter_effective_saturation(properties_array.effective_saturation, conc, geom_params, hydraulic_params);
% es(:, t_idx) = effective_saturation;
        if (mod(t_idx, progr) == 0)
            disp (t_idx / time_params.num_intervals * 100);
        end
    end
% mesh(es);
    plot(time_params.days_elapsed, smooth(leachate_flux, 1));
    
    toc;
    
    return
    
    function new_se = alter_effective_saturation(ini_se, add_conc, geom_params, hydraulic_params)
        column_volume_array = geom_params.dz .* geom_params.is_landfill_array .* geom_params.dx .* geom_params.dy;
        max_water_volume = (hydraulic_params.theta_s - hydraulic_params.theta_r) * column_volume_array;
        new_se = zeros(size(ini_se));
        idx = (max_water_volume ~= 0);
        new_se(idx) = ini_se(idx) + add_conc(idx) ./ max_water_volume(idx);
        % Make sure effective saturation doesn't exceed 1
        if ~isempty(find(new_se > 1, 1, 'first'))
            warning('Warning! Effective saturation exceeds 1');
        end
    end
end