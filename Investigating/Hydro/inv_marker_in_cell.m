function inv_marker_in_cell()
    %% Input data
    % Time
    time = struct();
    time.dt = 1;
    time.max_t = 1000;
    time.t = 0:time.dt:time.max_t;

    % Domain geometry
    geometry = struct();
    geometry.z_bot = -10;
    geometry.z_ref = -9;
    geometry.dz = 0.1;
    geometry.nz = (0 - geometry.z_bot) / geometry.dz;
    geometry.h_amb = geometry.z_ref;
    z = -geometry.dz * (0.5:geometry.nz)';
    h_ini =  geometry.h_amb - z;

    % Hydraulic properties
    hydro = struct();
    hydro.k_sat = 1;
    hydro.alpha   = 2;
    hydro.theta_r = 0.15;
    hydro.theta_s = 0.5;
    hydro.lambda  = 0.4;
    hydro.n       = hydro.lambda + 1;
    hydro.m       = hydro.lambda ./ hydro.n;

    % Markers
    num_markers = 1e+3;                     % initial number of markers
    dv = 1e-3;                              % volume of a single marker
    markers = nan(num_markers, 1);          % array of markers' coordinates
    
    %% Simulated data
    h = zeros(numel(time.t), geometry.nz);
    theta = zeros(numel(time.t), geometry.nz);
    
    %% Initial state
    [k_ini, theta_ini, ~] = van_genuchten(h_ini, hydro);
    h(1, :) = h_ini;
    theta(1, :) = theta_ini;
    
    %% Simulation
    v_in = -0.01;                           % volume of water in per 1m^2
    n_markers = v_in / dv;
    
    z_mark = zeros(n_markers, 1);           % position of markers
    
    u = v_in ./ theta_ini;                   % velocity field
    
    out = mic(z_mark, u, geometry);
    
    return
    
end