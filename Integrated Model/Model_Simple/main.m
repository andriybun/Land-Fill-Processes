function main()

%   Integrated modelling framework for landfills
%
%   Features:

    clc;
    addpath('../Common/')
    
    %% Input parameters:
    
    % Generate precipitation data (specific, independent of area):
    file_name = '../Data/precipitation_daily_KNMI_20110908.txt';
    [precipitation_intensity_time_vector, time_params, start_date] = read_precipitation_data_csv(file_name);
    
    % Fluid velocity parameters
    expected_fluid_velocity_mps = 1e-7;                     % m / s
    expected_fluid_velocity = ...
        expected_fluid_velocity_mps * time_params.time_discretization;  % m / {time step}
    variance = 3 * expected_fluid_velocity;                 % 95% confidence interval width

    %% Preparing data for simulation:
    
    % Get spatial characteristics of a landfill:
    spatial = define_geometry();
    
    % Determine probability distribution parameters corresponding to
    % defined inputs:
    [mu, sigma] = compute_lognormal_parameters(expected_fluid_velocity, variance);
    
    % Generate random flow velocities array:
    velocity_array = generate_velocities(spatial, spatial.column_height_array, mu, sigma);
    
%%% TODO: dynamically model residence times
    % Compute array of residual times for columns:
    velocity_non_zero_array = (velocity_array ~= 0);
    residual_time_array = zeros(size(velocity_array));
    residual_time_array(velocity_non_zero_array == 1) = spatial.column_height_array(velocity_non_zero_array == 1) ...
        ./ velocity_array(velocity_non_zero_array == 1);
    
%%% TODO: try different variants (floor, ceil, round)
    residual_time_array = ceil(residual_time_array / time_params.time_discretization);	% Discretized residual time in number of time intervals

    % Generate biogeochemical properties:
    properties_array = generate_biogeochemical_properties_2d(spatial);
    
    % Net amounts of precipitation from rainfall events entering each
    % column / pathway:
    precipitation_in_time_vector = precipitation_intensity_time_vector * (spatial.dx * spatial.dy);
    
    % Array specifying amounts of leachate leaving the landfill
	leachate_out_array = zeros(spatial.xn, spatial.yn, time_params.num_intervals);
    leachate_out = zeros(1, time_params.num_intervals);
    
    %% Main loop over time
    for t = 1:time_params.num_intervals
        if (precipitation_in_time_vector(t) > 0)
            leachate_out_array = set_out_flow(spatial, t, time_params.num_intervals, ...
                precipitation_in_time_vector(t), leachate_out_array, residual_time_array);
        end

        leachate_out(t) = sum(sum(leachate_out_array(:, :, t)));
        
        % Display progress
        if (mod(t, 1000) == 0)
            disp (t / time_params.num_intervals * 100);
        end
    end
    
    plot(leachate_out);
    
    data = struct();
    data.t = time_params.time_discretization : time_params.time_discretization : time_params.num_intervals * time_params.time_discretization;
    data.in = precipitation_in_time_vector;
    data.out = leachate_out;
    
    save('data.mat', '-struct', 'data');
    
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