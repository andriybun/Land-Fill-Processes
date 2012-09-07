function execute_richards_1d_loop()

    %%%%%%%%%%%
    %%%
    %%% first time step is void
    %%%
	%%%%%%%%%%%

    %% Loop over water table:
    data = struct();
    
    k_sat = 1;
    k_threshold = k_sat + 0.5;
    
    wt_min = 0.1;
    wt_max = 20;
    wt_nsteps = 20;
    data.water_table_elevation_vector = -10 .^ (linspace(log10(wt_min), log10(wt_max), wt_nsteps));
    plot(data.water_table_elevation_vector);

%     data.water_table_elevation_vector = [-0.1:-0.4:-19.9];
    num_time_steps = 100;
    data.t_range = zeros(numel(data.water_table_elevation_vector), num_time_steps);
    
    data.out_flux = zeros(numel(data.water_table_elevation_vector), num_time_steps);
    data.saturation_effective_avg = zeros(numel(data.water_table_elevation_vector), 1);

    if k_threshold >= 1
        % Matrix domain:
        t_max = 1e+3 * [ ...
            0.000060599876945
            0.000080799835927
            0.000111099774399
            0.000166649661599
            0.000249974492399
            0.000403999179634
            0.000719623538724
            0.001464497026174
            0.003226943447329
            0.007383085007817
            0.016180167144355
            0.032875433242743
            0.063023872022954
            0.114869591744400
            0.201423890986180
            0.344449700556220
            0.574941332517090
            0.939904091419243
            1.527594123048668
            2.450000000000000
            ];
    else
        % Channel domain:
        t_max = [ ...
            0.001391181748049
            0.002047819533136
            0.002785145859589
            0.004380058445863
            0.006730769160708
            0.010997549344592
            0.020264880796635
            0.041240810042243
            0.091859730823874
            0.212430361412549
            0.460592453145092
            0.945910843225541
            1.813359035860935
            3.305093854904683
            5.795483854515127
            9.910704576762988
            16.366557674958486
            27.043460236733655
            43.485232152900565
            70.492806856443238
            ];
    end;
    t_max = t_max / k_sat;
    
    for i = 1:numel(data.water_table_elevation_vector)
        max_t = t_max(i);
        
        data.t_range(i, :) = linspace(0, max_t, num_time_steps);
        dt = data.t_range(i, 2);
        [data.out_flux(i, :), data.saturation_effective_avg(i), t_max(i), vg_par, domain_name] = ...
            run_richards_fdm_implicit(data.water_table_elevation_vector(i), k_sat, data.t_range(i, :), dt, k_threshold);
        plot(data.t_range(i, :), data.out_flux(i, :));
    end

    % t_max .* data.water_table_elevation_vector(i)
    % [data.saturation_effective_avg t_max]
    % plot(data.saturation_effective_avg, t_max);
    
    data.t_range = data.t_range;
    data.van_genuchten_params = vg_par;
    save(sprintf('../Common/data_wt_loop_%s_domain.mat', domain_name), '-struct', 'data');

end