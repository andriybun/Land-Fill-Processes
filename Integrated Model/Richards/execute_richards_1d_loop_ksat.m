function execute_richards_1d_loop_ksat()

    %%%%%%%%%%%
    %%%
    %%% first time step is void
    %%%
	%%%%%%%%%%%

    % Loop over hydraulic conductivity:
    data = struct();
    
    time_factor = 1;
    
    data.water_table_elevation_vector = -2;
    
    k_sat_min = 0.1;
    k_sat_max = 300;
    k_sat_nsteps = 10;
    data.k_sat_vector = 10 .^ (linspace(log10(k_sat_min), log10(k_sat_max), k_sat_nsteps));
%     plot(data.k_sat_vector);
    
    num_time_steps = 100;
    data.t_range = zeros(numel(data.water_table_elevation_vector), num_time_steps);
    
    data.out_flux = zeros(numel(data.water_table_elevation_vector), num_time_steps);
    data.saturation_effective_avg = zeros(numel(data.water_table_elevation_vector), 1);

    t_max = zeros(size(data.k_sat_vector));
    
    for i = 1:numel(data.k_sat_vector)
        k_sat = data.k_sat_vector(i);

        max_t = 0.8 / k_sat;
        
        data.t_range(i, :) = linspace(0, max_t, num_time_steps);
        dt = data.t_range(i, 2);
        [data.out_flux(i, :), data.saturation_effective_avg(i), t_max(i), vg_par, domain_name] = ...
            run_richards_fdm_implicit(data.water_table_elevation_vector(1), k_sat, data.t_range(i, :), dt, k_sat_max + 1);
        plot(data.t_range(i, :), data.out_flux(i, :));
    end
    
    % t_max .* data.k_sat_vector
    
    data.t_range = data.t_range / time_factor;
    data.van_genuchten_params = vg_par;
    save(sprintf('../Common/data_ksat_loop_%s_domain.mat', domain_name), '-struct', 'data');
end