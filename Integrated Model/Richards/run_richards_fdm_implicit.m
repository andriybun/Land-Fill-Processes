% Script controlling Richards Equation
% using implicit finite difference solution
% Features:
%   * applying bottom pressure;
%   * applying flux on top;
%   * possibility to accumulate water above the system and flush it when
%     this is possible.

function [out_flux, saturation_effective_avg, t_max, vg_par, domain_name] = run_richards_fdm_implicit(water_table_elevation, k_sat, t_range, dt)
    %% Inputs

    tic;
    
    EPSILON = 1e-15;
    
    % Column size
    z_top = 0;
    z_bot = -1;
    
    % In flow into the system defined as a function below

    % Nodes/internodes
    nn = 20;
    zin = linspace(z_top, z_bot, nn+1)';
    dzin = zin(2:end) - zin(1:end-1);
    nin = nn + 1;
    zn(:, 1) = (zin(1:nin-1, 1) + zin(2:nin, 1)) / 2;
    dzn = zn(2:end) - zn(1:end-1);

    spatial.zn = zn;
    spatial.nn = nn;
    spatial.dzn = dzn;
    spatial.dzin = dzin;

    if nargin == 0
        % Position of water table, compare to the top of column
        water_table_elevation = -9.95;
    end
    
    % Water pressure
    hw = water_table_elevation - zn;
    
    % Water pressure below the system
    h_out = hw(end) - spatial.dzn(end);
%     h_out = -1.5-zn(end);
    
    % Time variables
    if nargin == 0
        t_range = [0, 20]; % days
        dt = 5e-2;
        t_range = t_range(1):dt:t_range(2);
    end
    num_time_steps = numel(t_range);

    if nargin == 0
        vg_par.k_sat = 0.1;   % m/day
    else
        vg_par.k_sat = k_sat;
    end
    
    % Van Genuchten parameters
    if vg_par.k_sat < 10
        % Matrix domain
        domain_name = 'matrix';
        vg_par.alpha   = 2;
        vg_par.theta_r = 0.15;
        vg_par.theta_s = 0.5;
        vg_par.lambda  = 0.4;
        vg_par.n       = vg_par.lambda + 1;
        vg_par.m       = vg_par.lambda ./ vg_par.n;
    else
        % Channel domain
        domain_name = 'channel';
        vg_par.alpha   = 2;
        vg_par.theta_r = 0;
        vg_par.theta_s = 0.01;
        vg_par.lambda  = 0.4;
        vg_par.n       = vg_par.lambda + 1;
        vg_par.m       = vg_par.lambda ./ vg_par.n;
    end
    
    %% Initializing variables
    clc();
    close all;

    % Accumulate water on top of the system and create pressure above:
    accumulate_water = 0;

    % Amount of water cumulated on top of column, which didn't enter the
    % system yet:
    hw_above = 0;
    dhw_above = 0;
    
    % Arrays for storing solution:
    theta = zeros(nn, num_time_steps);
    h = zeros(nn, num_time_steps);
    in_flux_cumulative = zeros(1, num_time_steps);
    in_precipitation_cumulative = zeros(1, num_time_steps);
    out_flux = zeros(1, num_time_steps);
    breakthrough = zeros(1, num_time_steps);

    % Initial state
    [k_n, theta(:, 1)] = van_genuchten(hw, vg_par);
    h(:, 1) = hw(1:nn);

    % Total relative moisture content
    saturation_effective_avg = mean((theta(:, 1) - vg_par.theta_r) / (vg_par.theta_s - vg_par.theta_r));
%     saturation_effective_avg = sum((theta(:, 1) - vg_par.theta_r) .* (-spatial.dzin)) / ((vg_par.theta_s - vg_par.theta_r) * (z_top - z_bot));
    
    %% Main time-loop
    for t_idx = 2:num_time_steps
        % Current time in seconds
        t = t_range(t_idx);
        dt = t_range(t_idx) - t_range(t_idx-1);
        in_flx = -1e-5 * vg_par.k_sat * (t_idx == 2);%-in_flux(t_idx);
        
        % Solve system of equations
        theta_prev = theta(:, t_idx-1);
        options = optimset('Display', 'off', 'TolFun', EPSILON);  % Turn off display % 
        % 'Algorithm', 'trust-region-reflective'
        % 'Algorithm', 'trust-region-reflective'
        % 'Algorithm', 'levenberg-marquardt'
        hw_new = fsolve(@(hw) net_flux(hw, t, dt, theta_prev, k_n, in_flx, hw_above, spatial, vg_par, h_out, accumulate_water, EPSILON), hw, options);
        [~, q_in, dhw_above] = net_flux(hw_new, t, dt, theta_prev, k_n, in_flx, hw_above, spatial, vg_par, h_out, accumulate_water, EPSILON);
        
        if accumulate_water
            hw_above = max(0, hw_above + dhw_above * dt);
        end

        hw = hw_new(1:spatial.nn);
        
        % Compute new value for hydraulic head
        [k_n, theta(:, t_idx)] = van_genuchten(hw, vg_par);
        
        % Save results
        h(:, t_idx) = hw;
        in_flux_cumulative(t_idx) = in_flux_cumulative(t_idx-1) - q_in(1) * dt;
        in_precipitation_cumulative(t_idx) = in_precipitation_cumulative(t_idx-1) - in_flx * dt;
        out_flux(t_idx) = q_in(spatial.nn+1) * dt;
        breakthrough(t_idx) = breakthrough(t_idx-1) - out_flux(t_idx);

        % Output
        fprintf('%5.2f\n', t_idx / num_time_steps * 100);
    end

    %% Displaying plots
    close all;
%     figure(1); hold on;
%     [x, y] = meshgrid(spatial.zn, t_range);
%     mesh(x, y, theta');
%     figure(2); hold on;
%     plot(t_range, breakthrough);
%     figure(2); hold on;
%     plot(t_range, in_flux_cumulative, 'r');
%     figure(2); hold on;
%     plot(t_range, in_precipitation_cumulative, 'g');
%     hold off;
%     figure(3); hold on;
%     plot(t_range, out_flux);
    
    t_max = t_range(find(breakthrough' / in_flux_cumulative(end) > 0.99, 1, 'first'));
    if isempty(t_max)
        warning(sprintf('%f\n', water_table_elevation));
        t_max = t_range(end);
    end

    toc;
    
end