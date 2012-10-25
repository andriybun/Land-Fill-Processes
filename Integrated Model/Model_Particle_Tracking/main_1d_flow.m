function main_1d_flow()
    clc;
	tic;
    addpath('../Common/');

    % Time params
    start_date = struct();
    start_date.year = 2012;
    start_date.month = 1;
    start_date.day = 1;
    time_params.max_days = 700;                                                            % number of simulation days
    time_params.intervals_per_day = 1;
    time_params.time_discretization = 1 / time_params.intervals_per_day;                 % in days
    time_params.num_intervals = time_params.max_days * time_params.intervals_per_day;    % in {time step}
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
    precipitation_intensity_time_vector = precipitation_intensity_time_vector * 1e-2 * time_params.time_discretization;
%     %% Input pulse:
%     precipitation_intensity_time_vector = zeros(size(time_params.days_elapsed));
%     precipitation_intensity_time_vector(:) = 0;
%     precipitation_intensity_time_vector(1) = 1e-4 * time_params.time_discretization;
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
%     plot(time_params.days_elapsed, leachate_flux(:, end), 'Color', [0.6, 0.4, 0], 'LineWidth', 1.7);
    plot(time_params.days_elapsed, leachate_flux(:, 1:9:10), 'LineWidth', 1.7);
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
    
    function [leachate_intercell_array, properties_array] = transport_lognormal(leachate_intercell_array, t, properties_array, scale, ...
            geometry_params, hydraulic_params, time_params, lognrnd_param_definer)
%         [mu, sigma] = lognrnd_param_definer.get_params(...
%             hydraulic_params.k_sat ./ geometry_params.column_height_array, properties_array.effective_saturation);
        [mu, sigma] = lognrnd_param_definer.get_params(...
            hydraulic_params.k_sat, properties_array.effective_saturation); 
        mu = 3.5075 * log(geometry_params.column_height_array) + mu;
        
        pow = 1.7;
        
        steps = power(linspace(0, 1, geometry_params.zn + 1), pow)';
        mean_orig = exp(mu + sigma * sigma / 2);
        mu_steps = log(mean_orig * steps(2:end)) - sigma * sigma / 2;
        
        sz = size(geometry_params.is_landfill_array);
        if numel(sz) == 2
            sz = [sz, 1];
        end
        
        for x_idx = 1:sz(2)
            for y_idx = 1:sz(3)
                z_idx = 1;
                while (z_idx <= geometry_params.zn)
                    z_idx_start = z_idx;
                    shift_flag = false;
                    while (z_idx <= geometry_params.zn) && (geometry_params.is_landfill_array(z_idx, x_idx, y_idx) == 0)
                        z_idx = z_idx + 1;
                        shift_flag = true;
                    end
                    if shift_flag
                        num_shift = geometry_params.zn - z_idx + 1;
                        mu_steps(z_idx:end) = mu_steps(z_idx_start:z_idx_start + num_shift - 1);
                        if z_idx_start == 1
                            mu_steps(z_idx_start:z_idx-1) = 0;
                        else
                            mu_steps(z_idx_start:z_idx-1) = mu_steps(z_idx_start-1);
                        end
                    end
                    z_idx = z_idx + 1;
                end
            end
        end

%         % Limiting influx        
%         se = properties_array.effective_saturation(1, :, :);
%         k = hydraulic_params.k_sat .* se .^ (1/2) .* (1 - (1 - se .^ (1 ./ hydraulic_params.m)) .^ hydraulic_params.m) .^ 2;
%         in_flux_lim = k * time_params.time_discretization;
%         scale = min(scale, in_flux_lim);
        
        if (t == 1)
            in_flux_cumulative(t) = scale;
        else
            in_flux_cumulative(t) = in_flux_cumulative(t-1) + scale;
        end
        
        t_vector = 0 : time_params.time_discretization : time_params.time_discretization * (num_intervals - t + 1);
        
        for cell_idx = 1:geometry_params.zn
            scale = scale * geometry_params.dx * geometry_params.dy;
            breakthrough_cum = scale * log_normal_cdf(t_vector, mu_steps(cell_idx), sigma);
            breakthrough = breakthrough_cum(2:end) - breakthrough_cum(1:end-1);
            leachate_intercell_array(t:end, cell_idx) = leachate_intercell_array(t:end, cell_idx) + breakthrough';
        end
        
        intercell_flux = scale - leachate_intercell_array(t, end);
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
        dt = t(2:end) - t(1:end-1);
        mn = sum(t(2:end) .* (out_flx(1:end-1) + out_flx(2:end)) ./ 2 .* dt);
        varn = sum((t(2:end) - mn).^2 .* (out_flx(1:end-1) + out_flx(2:end)) ./ 2 .* dt);
        [mux, sigmax] = get_mu_sigma(mn, varn);
    end

    function [mux, sigmax] = get_mu_sigma(evn, varn)
        sigmax = sqrt(log(varn ./ (evn .* evn) + 1));
        mux = log(evn) - sigmax .* sigmax ./ 2;
    end
end