function main_1d_flow_solute()
    clc;
    tic;
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
    
    % Flow params
    spatial_params = struct();
    spatial_params.dx = 1;
    spatial_params.dy = 1;
    spatial_params.dz = 1;                              % vertical dimensions of spatial_params.dz, forming column
    spatial_params.zn = 8;
    spatial_params.is_landfill_array = ones(spatial_params.zn, 1);
    spatial_params.is_landfill_array(5:end) = 0;
    
    file_name = '../Data/precipitation_daily_KNMI_20110908.txt';
    [precipitation_intensity_time_vector, time_params, start_date] = read_precipitation_data_csv(file_name);
%     precipitation_intensity_time_vector = zeros(1, num_intervals);
%     precipitation_intensity_time_vector(1) = 1e-4;
    num_intervals = time_params.num_intervals;

    % Determine probability distribution parameters corresponding to defined inputs:
    lognrnd_param_definer = log_normal_params('../Common/opt_params_wt_channel_domain.mat');
    
    % Fluid hydraulic parameters
    hydraulic_params = lognrnd_param_definer.hydraulic_params;
    hydraulic_params.k_sat_ref = hydraulic_params.k_sat;    % reference conductivity
    hydraulic_params.k_sat = 1e-3;                          % relative conductivity compared to reference conductivity
    hydraulic_params.d = 1e-1;                              % diffusion_coefficient

	% initial SE
    properties_array = generate_biogeochemical_properties_3d(spatial_params, hydraulic_params);

    % Initial solute to be flushed out
    l = spatial_params.is_landfill_array * spatial_params.dz;
    domain_type = 2;                                        % 1 - matrix, 2 - channel
    stc_ini = solute_transport_class(l, hydraulic_params.d, properties_array.solutes.solute_1_fraction, domain_type, 1);
    
    % Solute that is transported from other cells
    stc_trans = solute_transport_class(l, hydraulic_params.d, zeros(size(spatial_params.is_landfill_array)), domain_type, 2);
    stc_trans = repmat(stc_trans, [1, time_params.num_intervals]);
    
    % Result
    leachate_flux = zeros(numel(spatial_params.is_landfill_array), num_intervals);
    solute_out = zeros(numel(spatial_params.is_landfill_array), num_intervals);
    solute_mobile_out = zeros( numel(spatial_params.is_landfill_array), num_intervals);

    progr = floor(time_params.num_intervals / 10);
    
    for t_idx = 1:num_intervals
        [leachate_flux, properties_array] = transport_lognormal(leachate_flux, t_idx, properties_array, ...
            precipitation_intensity_time_vector(t_idx), spatial_params, hydraulic_params, time_params, lognrnd_param_definer);
        % Display progress
        if (mod(t_idx, progr) == 0)
            disp (t_idx / time_params.num_intervals * 100);
        end
    end

    plot(time_params.days_elapsed, leachate_flux(end, :), 'b');
    hold on;
    st = 1;
    plot(time_params.days_elapsed(st:end), -1e-1 * solute_out(4, st:end), 'r');
    hold off;

    toc;
    
    return 
    
    function [leachate_intercell_array, properties_array] = transport_lognormal(leachate_intercell_array, t_idx, ...
                            properties_array, scale, spatial_params, hydraulic_params, time_params, lognrnd_param_definer)
        [mu, sigma] = lognrnd_param_definer.get_params(...
            hydraulic_params.k_sat ./ spatial_params.dz, properties_array.effective_saturation);
        skip_cell_idx = (spatial_params.is_landfill_array == 0);
        mu(skip_cell_idx) = 0;
        sigma(skip_cell_idx) = Inf;
        
        t_vector = 0 : time_params.time_discretization : time_params.time_discretization * (num_intervals - t_idx + 1);
        
        scale = scale * spatial_params.dx * spatial_params.dy;
        breakthrough_cum = scale * log_normal_cdf(t_vector, mu(1), sigma(1));
        breakthrough = breakthrough_cum(2:end) - breakthrough_cum(1:end-1);
        leachate_intercell_array(1, t_idx:end) = leachate_intercell_array(1, t_idx:end) + breakthrough;
        for cell_idx = 2:spatial_params.zn
            breakthrough_cum = leachate_intercell_array(cell_idx - 1, t_idx) ...
                * log_normal_cdf(t_vector, mu(cell_idx), sigma(cell_idx));
            breakthrough = breakthrough_cum(2:end) - breakthrough_cum(1:end-1);
            leachate_intercell_array(cell_idx, t_idx:end) = leachate_intercell_array(cell_idx, t_idx:end) + breakthrough;
        end
        
        intercell_flux = vertcat(scale, leachate_intercell_array(1:end-1, t_idx)) - leachate_intercell_array(:, t_idx);
        properties_array.effective_saturation = alter_effective_saturation(properties_array.effective_saturation, ...
            intercell_flux, spatial_params, hydraulic_params);
        
        %% Solute transport module:
        % Average water velocity per cell - depends on conductivity:
        mean_water_velocity = exp(mu + sigma .^ 2 / 2);
        mean_water_velocity(mean_water_velocity == Inf) = 0;
        % Actual velocity is assumed to be a fraction of the average
        % velocity corresponding to the effective saturation:
        water_velocity = 1e-1 .* properties_array.effective_saturation .* mean_water_velocity;
        
        t = time_params.time_discretization * (t_idx - 1);
        if t_idx > 1
            % Flushing solute particles present from the beginning
            [stc_ini, solute_out(:, t_idx)] = stc_ini.flush(t, water_velocity);
            % Flushing mobile solute
            stc_trans(t_idx) = stc_trans(t_idx).set_c_ini(shift_array(-solute_out(:, t_idx) - solute_mobile_out(:, t_idx - 1)));
            % Loop over all mobile solute particles that were mobilized before
            for j = 1:t_idx
                [stc_trans(j), tmp_out] = stc_trans(j).flush(t, water_velocity);
                solute_mobile_out(:, t_idx) = solute_mobile_out(:, t_idx) + tmp_out;
            end
            
        else
            % Flushing solute particles present from the beginning
            [stc_ini, solute_out(:, t_idx)] = stc_ini.flush(t, water_velocity);
            % Flushing mobile solute
            stc_trans(t_idx) = stc_trans(t_idx).set_c_ini(shift_array(-solute_out(:, t_idx)));
            j = 1;
            [stc_trans(j), solute_mobile_out(:, t_idx)] = stc_trans(j).flush(t, water_velocity);
        end
    end
        
    function new_se = alter_effective_saturation(current_se, flx, spatial_params, hydraulic_params)
        column_volume_array = spatial_params.dz .* spatial_params.is_landfill_array .* spatial_params.dx .* spatial_params.dy;
        max_water_volume = (hydraulic_params.theta_s - hydraulic_params.theta_r) * column_volume_array;
        new_se = zeros(size(current_se));
        idx = (max_water_volume ~= 0);
        new_se(idx) = current_se(idx) + flx(idx) ./ max_water_volume(idx);
        % Make sure effective saturation doesn't exceed 1
        if ~isempty(find(new_se > 1, 1, 'first'))
            warning('Warning! Effective saturation exceeds 1');
        end
    end
    
    function c = shift_array(c_in)
        sz = size(c_in);
        if sz(1) == 1
            c = cat(2, 0, c_in(1:end-1));
        elseif sz(2) == 1
            c = cat(1, 0, c_in(1:end-1));
        else
            error('Error! Wrong size of array');
        end
    end
end