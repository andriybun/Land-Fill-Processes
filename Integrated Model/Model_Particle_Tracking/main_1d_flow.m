function main_1d_flow()
    clc;
	tic;
    addpath('../Common/');

    % Time params
    start_date = struct();
    start_date.year = 2012;
    start_date.month = 1;
    start_date.day = 1;
    time_params.max_days = 700;                                                             % number of simulation days
    time_params.intervals_per_day = 1;
    time_params.dt = 1 / time_params.intervals_per_day;                                     % in days
    time_params.num_intervals = time_params.max_days * time_params.intervals_per_day;       % in {time step}
    time_params.days_elapsed = (0 : 1: (time_params.num_intervals-1)) / time_params.intervals_per_day;

    % Geometry params
    geometry_params = struct();
    geometry_params.dx = 1;
    geometry_params.dy = 1;
    geometry_params.dz = 1;                              % vertical dimensions of geometry_params.dz, forming column
    geometry_params.zn = 10;
    geometry_params.is_landfill_array = ones(geometry_params.zn, 1);
    geometry_params.column_height_array = geometry_params.dz * squeeze(sum(geometry_params.is_landfill_array, 1));
    
    file_name = '../Data/precipitation_daily_KNMI_20110908.txt';
    [precipitation_intensity_time_vector, time_params, start_date] = read_precipitation_data_csv(file_name);
    precipitation_intensity_time_vector = precipitation_intensity_time_vector * 1e-2 * time_params.dt;
%     %% Input pulse:
%     precipitation_intensity_time_vector = zeros(size(time_params.days_elapsed));
%     precipitation_intensity_time_vector(:) = 0;
%     precipitation_intensity_time_vector(1) = 1e-4 * time_params.dt;
%     %% End input pulse
    num_intervals = time_params.num_intervals;
    
    % Determine probability distribution parameters corresponding to defined inputs:
    lognrnd_param_definer = log_normal_params('../Common/opt_params_wt_channel_domain.mat');
    
    % Fluid hydraulic parameters
    hydraulic_params = lognrnd_param_definer.hydraulic_params;
    hydraulic_params.k_sat_ref = hydraulic_params.k_sat;    % reference conductivity
    hydraulic_params.k_sat = 1e+4;                          % relative conductivity compared to reference conductivity
    hydraulic_params.d = 1;                                 % diffusion_coefficient

	% initial SE
    properties_array = generate_biogeochemical_properties_3d(geometry_params, hydraulic_params);
    
    % Result
    leachate_flux = zeros(num_intervals, numel(geometry_params.is_landfill_array));
    in_flux_cumulative = zeros(num_intervals, 1);

    progr = floor(time_params.num_intervals / 10);
    
    for t = 1:num_intervals
        [leachate_flux, properties_array] = transport_lognormal(leachate_flux, t, properties_array, ...
            precipitation_intensity_time_vector(t), geometry_params, hydraulic_params, time_params, lognrnd_param_definer);
        if (mod(t, progr) == 0)
            disp (t / time_params.num_intervals * 100);
        end
    end

    % First
    close all;
    figure('OuterPosition', [500, 600, 600, 400]);
    hold on;
    subplot(1, 1, 1);
    subplot('Position', [0.1 0.27 0.85 0.63]);
%     plot(time_params.days_elapsed, leachate_flux(:, end), 'Color', [0.6, 0.4, 0], 'LineWihydro_1d_2domainh', 1.7);
    plot(time_params.days_elapsed, leachate_flux(:, end), 'LineWidth', 1.7); % , 'Color', [0.8, 0, 0.8]
    xlabel('Days');
    ylabel('Out flux');
    annotation('textbox', [0.1, 0.05, 0.85, 0.075], 'String', ...
        sprintf('K_s_a_t = %3.3f', hydraulic_params.k_sat));
    hold off;
    % End first

%     %% Precipitation
%     x = 45;
%     figure('OuterPosition', [500, 600+x, 600, 400-x]);
%     subplot(1, 1, 1);
%     subplot('Position', [0.1 0.15 0.85 0.73]);
%     plot(time_params.days_elapsed, precipitation_intensity_time_vector, 'b');
%     xlabel('Days');
%     ylabel('Precipitation');
%     %% End precipitation
    
	toc;

    return 
    
    function [leachate_intercell_array, properties_array] = transport_lognormal(leachate_intercell_array, t_idx, properties_array, scale, ...
            geometry_params, hydraulic_params, time_params, lognrnd_param_definer)
        [mu, sigma] = lognrnd_param_definer.get_params(...
            hydraulic_params.k_sat, properties_array.effective_saturation); 
        mu = 3.5075 * log(geometry_params.column_height_array) + mu;

        log_normal_params = struct('mu', [], 'sigma', []);
        log_normal_params.mu = mu;
        log_normal_params.sigma = sigma;
        
        % Picked function that defines distribution of mean average flow over cells
        pow = 1.7;
        mean_steps_rel = power(linspace(0, 1, geometry_params.zn + 1), pow);

        % Calculate mu for intermediate points
        
        mu_steps = get_intermediate_mu(geometry_params, log_normal_params, mean_steps_rel);
        
%         % Limiting influx        
%         se = properties_array.effective_saturation(1, :, :);
%         k = hydraulic_params.k_sat .* se .^ (1/2) .* (1 - (1 - se .^ (1 ./ hydraulic_params.m)) .^ hydraulic_params.m) .^ 2;
%         in_flux_lim = k * time_params.dt;
%         scale = min(scale, in_flux_lim);
        scale = scale * geometry_params.dx * geometry_params.dy;
        
        % Calculating total flux that entered system
        if (t_idx == 1)
            in_flux_cumulative(t_idx) = scale;
        else
            in_flux_cumulative(t_idx) = in_flux_cumulative(t_idx-1) + scale;
        end
        
        % Flux simulation
        for cell_idx = 1:geometry_params.zn
            log_normal_params.mu = mu_steps(cell_idx);
            breakthrough = hydro_1d(scale, t_idx, geometry_params, time_params, log_normal_params);
            leachate_intercell_array(:, cell_idx) = leachate_intercell_array(:, cell_idx) + breakthrough';
        end
        
        intercell_flux = scale - leachate_intercell_array(t_idx, end);
        properties_array.effective_saturation = alter_effective_saturation(properties_array.effective_saturation, ...
            intercell_flux, geometry_params, hydraulic_params);
    end
        
    function new_se = alter_effective_saturation(current_se, flx, geometry_params, hydraulic_params)
        column_volume_array = geometry_params.column_height_array .* geometry_params.dx .* geometry_params.dy;
        max_water_volume = (hydraulic_params.theta_s - hydraulic_params.theta_r) * column_volume_array;
        new_se = zeros(size(current_se));
        idx = (max_water_volume ~= 0);
        new_se(idx) = current_se(idx) + flx(idx)' ./ max_water_volume(idx);
        % Make sure effective saturation doesn't exceed 1
        if ~isempty(find(new_se > 1, 1, 'first'))
            warning('Warning! Effective saturation exceeds 1');
        end
    end
   
    function [evx, varx] = get_moments(mux, sigmax)
        % Expected value
        evx = exp(mux + sigmax .* sigmax ./ 2);
        % Variance
        varx = (exp(sigmax .* sigmax) - 1) .* exp(2 .* mux + sigmax .* sigmax);
    end

    function [mux, sigmax] = pick_params(t, out_flx)
        hydro_1d_2domain = t(2:end) - t(1:end-1);
        mn = sum(t(2:end) .* (out_flx(1:end-1) + out_flx(2:end)) ./ 2 .* hydro_1d_2domain);
        varn = sum((t(2:end) - mn).^2 .* (out_flx(1:end-1) + out_flx(2:end)) ./ 2 .* hydro_1d_2domain);
        [mux, sigmax] = get_mu_sigma(mn, varn);
    end

    function [mux, sigmax] = get_mu_sigma(evn, varn)
        sigmax = sqrt(log(varn ./ (evn .* evn) + 1));
        mux = log(evn) - sigmax .* sigmax ./ 2;
    end
end