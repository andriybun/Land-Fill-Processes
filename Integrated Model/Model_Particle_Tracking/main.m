function main()

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
    expected_fluid_velocity_mps = 1e-7;                     % m / s
    expected_fluid_velocity = ...
        expected_fluid_velocity_mps * time_params.time_discretization;  % m / {time step}
    variance = 3 * expected_fluid_velocity;                 % 95% confidence interval width

    %% Preparing data for simulation:
    
    % Get spatial characteristics of a landfill:
    [spatial, is_landfill_array, column_height_array] = define_geometry();
    
    % Determine probability distribution parameters corresponding to
    % defined inputs:
    [mu, sigma] = compute_lognormal_parameters(expected_fluid_velocity, variance);
    
    % Generate random flow velocities array:
    velocity_array = generate_velocities(spatial, column_height_array, mu, sigma);
    
    % Compute array of residual times for columns:
    velocity_non_zero_array = (velocity_array ~= 0);
    residual_time_array = zeros(size(velocity_array));
    residual_time_array(velocity_non_zero_array == 1) = column_height_array(velocity_non_zero_array == 1) ./ velocity_array(velocity_non_zero_array == 1);
    
%%% TODO: try different variants (floor, ceil, round)
    % Array of residence times per column
    residual_time_array = ceil(residual_time_array ./ time_params.time_discretization);	% Discretized residual time in number of time intervals
    % Array of numbers of cells passed in one time interval
    num_cells_passed = velocity_array .* time_params.time_discretization ./ spatial.dz;

    % Generate biogeochemical properties:
    properties_array = generate_biogeochemical_properties(spatial, column_height_array);
    
    % Generate precipitation data (specific, independent of area):
    precipitation_intensity_time_vector = generate_precipitation_data(start_date, time_params);
    % Net amounts of precipitation from rainfall events entering each
    % column / pathway:
    precipitation_in_time_vector = precipitation_intensity_time_vector * (spatial.dx * spatial.dy);
    
    % Array specifying amounts of leachate leaving the landfill
	leachate_out_array = zeros(spatial.xn, spatial.yn, num_intervals);
    leachate_out = zeros(1, num_intervals);
    
    %% Main loop over time
    for t = 1:num_intervals
        if (precipitation_in_time_vector(t) > 0)
            leachate_out_array = set_out_flow(spatial, t, num_intervals, precipitation_in_time_vector(t), leachate_out_array, residual_time_array);
        end

        leachate_out(t) = sum(sum(leachate_out_array(:, :, t)));
        
        % Display progress
        if (mod(t, 1000) == 0)
            disp (t / num_intervals * 100);
        end
    end
    
    plot(leachate_out);
    
    return
    
    function leachate_out_array = set_out_flow(spatial, t, num_intervals, in_flow, leachate_out_array, residual_time_array)
        time_out = t + residual_time_array;
        test_idx = (time_out <= num_intervals) .* (residual_time_array > 0);
        
        for rw = 1:spatial.xn;
            for cl = 1:spatial.yn
                if (residual_time_array(rw, cl) > 0) && (test_idx(rw, cl))
                    leachate_out_array(rw, cl, time_out(rw, cl)) = in_flow;
                end
            end
%             idx = (residual_time_array(rw, :) > 0) .* (test_idx(rw, :));
%             leachate_out_array(rw, idx, time_out(rw, idx)) = in_flow;
%             squeeze(leachate_out_array(rw, :, :))
        end
    end
end