function execute_richards_1d_loop()

    %%%%%%%%%%%
    %%%
    %%% first time step is void
    %%%
	%%%%%%%%%%%

    %% Loop over water table:
    data = struct();
    
    k_sat = 1;
    k_threshold = k_sat - 0.5;
    
    wt_min = 0.1;
    wt_max = 40;
    wt_nsteps = 25;
    data.water_table_elevation_vector = -10 .^ (linspace(log10(wt_min), log10(wt_max), wt_nsteps));
    data.water_table_elevation_vector(9) = -1;
    plot(data.water_table_elevation_vector);

%     data.water_table_elevation_vector = [-0.1:-0.4:-19.9];
    num_time_steps = 100;
    data.t_range = zeros(numel(data.water_table_elevation_vector), num_time_steps);
    
    data.out_flux = zeros(numel(data.water_table_elevation_vector), num_time_steps);
    data.saturation_effective_avg = zeros(numel(data.water_table_elevation_vector), 1);

    if k_threshold >= 1
        % Matrix domain:
        t_max = 1.5e+3 * [ ...
            0.000047263081254
            0.000066249105347
            0.000093314288629
            0.000131299733381
            0.000192460720298
            0.000290849730792
            0.000470057126305
            0.000840655062672
            0.003972071706531
            0.003363480722365
            0.007057465007868
            0.014168536885761
            0.027161803430768
            0.049506086677721
            0.085904689758416
            0.145075845719424
            0.236523156076819
            0.377457657772282
            0.592271979143449
            0.926155494337312
            1.440495867768596
            2.184879093970003
            3.309891847770635
            5.002632384450566
            7.616187123762882
            ];
    else
        % Channel domain:
        t_max = 1.5e+2 * [ ...
            0.000013630770663
            0.000018838017458
            0.000026726147137
            0.000037226698147
            0.000055069929497
            0.000083314767762
            0.000134375533412
            0.000238330729598
            0.001133877983250
            0.000957139553104
            0.002009012077538
            0.004076652792353
            0.007733742240484
            0.014095798707787
            0.024456811776629
            0.041307242318209
            0.067329824043906
            0.107472990125577
            0.170373728667624
            0.263673737308153
            0.407878521169373
            0.624859310002755
            0.949196447318161
            1.436272174397226
            2.169737931765925
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