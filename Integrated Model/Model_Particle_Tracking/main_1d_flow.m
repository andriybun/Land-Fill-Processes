function main_1d_flow()
    clc;
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
    
    % Flow params
    spatial_params = struct();
    spatial_params.dx = 1;
    spatial_params.dy = 1;
    spatial_params.dz = 1;                              % vertical dimensions of spatial_params.dz, forming column
    spatial_params.zn = 25;
    spatial_params.is_landfill_array = ones(spatial_params.zn, 1);
%     spatial_params.is_landfill_array(5:end) = 0;
    
    file_name = '../Data/precipitation_daily_KNMI_20110908.txt';
%     [precipitation_intensity_time_vector, time_params, start_date] = read_precipitation_data_csv(file_name);
%     precipitation_intensity_time_vector = precipitation_intensity_time_vector * 1e-2 * time_params.time_discretization;
    precipitation_intensity_time_vector = zeros(size(time_params.days_elapsed));
    precipitation_intensity_time_vector(:) = 0;
    precipitation_intensity_time_vector(1) = 1e-4 * time_params.time_discretization;
    num_intervals = time_params.num_intervals;
    
    % Determine probability distribution parameters corresponding to defined inputs:
    lognrnd_param_definer = log_normal_params('../Common/opt_params_wt_channel_domain.mat');
    
    % Fluid hydraulic parameters
    hydraulic_params = lognrnd_param_definer.hydraulic_params;
    hydraulic_params.k_sat_ref = hydraulic_params.k_sat;    % reference conductivity
    hydraulic_params.k_sat = 1e-0;                          % relative conductivity compared to reference conductivity
    hydraulic_params.d = 1;                                 % diffusion_coefficient

	% initial SE
    properties_array = generate_biogeochemical_properties_3d(spatial_params, hydraulic_params);
    
    % Result
    leachate_flux = zeros(num_intervals, numel(spatial_params.is_landfill_array));
    in_flux_cumulative = zeros(num_intervals, 1);

    progr = floor(time_params.num_intervals / 10);
    
    for t = 1:num_intervals
        [leachate_flux, properties_array] = transport_lognormal(leachate_flux, t, properties_array, ...
            precipitation_intensity_time_vector(t), spatial_params, hydraulic_params, time_params, lognrnd_param_definer);
        % Display progress
        if (mod(t, progr) == 0)
            disp (t / time_params.num_intervals * 100);
        end
    end

%     % First
%     close all;
%     figure('OuterPosition', [500, 600, 600, 400]);
%     hold on;
%     subplot(1, 1, 1);
%     subplot('Position', [0.1 0.27 0.85 0.63]);
%     plot(time_params.days_elapsed, leachate_flux(:, end), 'b', 'LineWidth', 1.7); %  / time_params.time_discretization
%     xlabel('Days');
%     ylabel('Out flux');
%     annotation('textbox', [0.1, 0.05, 0.85, 0.075], 'String', ...
%         sprintf('K_s_a_t = %3.3f', hydraulic_params.k_sat));
%     hold off;
%     % End first

    %% Second
    hold on;
    plot(time_params.days_elapsed, leachate_flux(:, end), 'Color', [0, 0.8, 0], 'LineWidth', 1.3);
    hold off;
    
    [mu, sigma] = pick_params(time_params.days_elapsed', leachate_flux(:, end) / sum(leachate_flux(:, end) * time_params.time_discretization));
    disp([spatial_params.zn, mu, sigma]);
    
    %% End second
    
%     %% Precipitation
%     x = 45;
%     figure('OuterPosition', [500, 600+x, 600, 400-x]);
%     subplot(1, 1, 1);
%     subplot('Position', [0.1 0.15 0.85 0.73]);
%     plot(time_params.days_elapsed, precipitation_intensity_time_vector, 'b');
%     xlabel('Days');
%     ylabel('Precipitation');
%     %% End precipitation
    
    return 
    
    function [leachate_intercell_array, properties_array] = transport_lognormal(leachate_intercell_array, t, properties_array, scale, ...
            spatial_params, hydraulic_params, time_params, lognrnd_param_definer)
        [mu, sigma] = lognrnd_param_definer.get_params(...
            hydraulic_params.k_sat ./ spatial_params.dz, properties_array.effective_saturation);
%         [ev, vr] = get_moments(mu, sigma);
%         [mu_cum, sigma_cum] = get_mu_sigma(cumsum(flipud(ev)), cumsum(flipud(vr)));
%         [cumsum(spatial_params.is_landfill_array) * spatial_params.dz, mu_cum, sigma_cum];
        
        se = properties_array.effective_saturation(1, :, :);
        k = hydraulic_params.k_sat .* se .^ (1/2) .* (1 - (1 - se .^ (1 ./ hydraulic_params.m)) .^ hydraulic_params.m) .^ 2;
        
        % Limiting influx
        in_flux_lim = k * time_params.time_discretization;
        scale = min(scale, in_flux_lim);
        if (t == 1)
            in_flux_cumulative(t) = scale;
        else
            in_flux_cumulative(t) = in_flux_cumulative(t-1) + scale;
        end
        
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