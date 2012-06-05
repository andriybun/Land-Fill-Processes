function main_transfer_function()

%   Integrated modelling framework for landfills
%
%   Features:
%     - mass exchange between pathways;

    clc;
    addpath('../Common/')

    %% Input parameters:
    
    % Time
    start_date = struct();
    start_date.month = 1;
    start_date.day = 1;
    time_params.max_days = 30; % 365;                                                           % number of simulation days
    time_params.time_discretization = 3600;                                                     % in seconds
    time_params.intervals_per_day = 24 * 3600 / time_params.time_discretization;
    num_intervals = time_params.max_days * time_params.intervals_per_day;                       % in {time step}
    time_params.num_intervals = num_intervals;
    
    % Fluid velocity parameters
    k_sat = 2e-5;                                           % m / s
%     expected_fluid_velocity = ...
%         expected_fluid_velocity_mps * time_discretization;  % m / {time step}
%     variance = 3 * expected_fluid_velocity;                 % 95% confidence interval width

    %% Preparing data for simulation:
    
    % Get spatial characteristics of a landfill:
    [spatial, is_landfill_array, column_height_array] = define_geometry();
    spatial.is_landfill_array = is_landfill_array;
    spatial.column_height_array = column_height_array;
    
    % Determine probability distribution parameters corresponding to defined inputs:
    param_definer = log_normal_params();

    % Generate biogeochemical properties:
    properties_array = generate_biogeochemical_properties(spatial, column_height_array);
    
    % Generate precipitation data (specific, independent of area):
    rand('seed', 1);
    precipitation_intensity_time_vector = generate_precipitation_data(start_date, time_params);
    % Net amounts of precipitation from rainfall events entering each column / pathway:
    precipitation_in_time_vector = precipitation_intensity_time_vector * (spatial.dx * spatial.xn * spatial.dy * spatial.yn);
    
    % Array specifying amounts of leachate leaving the landfill
	leachate_out_array = zeros(num_intervals, spatial.xn, spatial.yn);
    leachate_out = zeros(1, num_intervals);
    
    %% Main loop over time
    for t = 1:num_intervals
        if (precipitation_in_time_vector(t) > 0)
            leachate_out_array(t:end, :, :) = leachate_out_array(t:end, :, :) + ...
                transport_lognormal(t, precipitation_in_time_vector(t), spatial, time_params, param_definer);
        end

        % Display progress
        if (mod(t, 1000) == 0)
            disp (t / num_intervals * 100);
        end
    end

    leachate_out = squeeze(sum(sum(leachate_out_array, 3), 2));
    plot(leachate_out);
    
    return
    
    function breakthrough = transport_lognormal(t, scale, spatial, time_params, param_definer)
        num_intervals = time_params.num_intervals;
        time_discretization = time_params.time_discretization;
        
        t_vector = 0 : time_discretization : time_discretization * (num_intervals - t);
        z = spatial.dz * spatial.zn;
        %% TODO: change 0.3 to dynamic value individual for each cell
        [mu, sigma] = param_definer.get_params(spatial.column_height_array * k_sat, 0.3);
        breakthrough = zeros(num_intervals - t + 1, size(mu, 1), size(mu, 2));
        % Leave breakthrough(1, :, :) = 0;
        for t_idx = 2:num_intervals - t + 1
            breakthrough(t_idx, :, :) = scale * exp(-(log(t_vector(t_idx)) - mu).^2 ./ ...
                (2 .* sigma  .* sigma)) ./ (sqrt(2 * pi) .* sigma .* t_vector(t_idx));
        end
    end
end