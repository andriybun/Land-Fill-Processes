function execute_richards_1d_loop_vgpar()

    %%%%%%%%%%%
    %%%
    %%% first time step is void
    %%%
	%%%%%%%%%%%

    % Loop over hydraulic conductivity:
    data = struct();
    
    time_factor = 1;
    
    data.water_table_elevation_vector = -2;
    
    k_sat = 1;

    shock_vector = roundn(exp(linspace(log(0.2), log(30), 15)), -6);
%     0:0.5:10;
    
    num_time_steps = 100;
    data.t_range = zeros(numel(shock_vector), num_time_steps);
    
    data.out_flux = zeros(numel(shock_vector), num_time_steps);
    data.saturation_effective_avg = zeros(numel(shock_vector), 1);

    t_max = zeros(size(shock_vector));
    data.length_vector = zeros(size(shock_vector));
    
    t_max = 1.0e+03 * [ ...
        0.000123420056262
        0.000246366979778
        0.000490010891129
        0.000935601807402
        0.001826491754630
        0.003443169940172
        0.006260067638514
        0.011211809646961
        0.020727314254609
        0.042160110655918
        0.100125389566483
        0.265323801939863
        0.722856586732834
        2.042232663057813
        5.910032549080127
        ];

    for i = 1:numel(shock_vector)
        max_t = t_max(i);
        data.t_range(i, :) = linspace(0, max_t, num_time_steps);
        dt = data.t_range(i, 2);
        [data.out_flux(i, :), data.saturation_effective_avg(i), t_max(i), vg_par, domain_name] = ...
            run_richards_fdm_implicit(data.water_table_elevation_vector(1), k_sat, data.t_range(i, :), dt, 0.5, shock_vector(i));
        data.length_vector(i) = shock_vector(i);
        plot(data.t_range(i, :), data.out_flux(i, :));
    end
    
    % t_max .* data.k_sat_vector
    
    data.t_range = data.t_range / time_factor;
    data.van_genuchten_params = vg_par;
    save(sprintf('../Common/data_length_loop_%s_domain.mat', domain_name), '-struct', 'data');
end