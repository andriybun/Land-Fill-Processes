function execute_richards_1d_loop_vgpar()

    %%%%%%%%%%%
    %%%
    %%% first time step is void
    %%%
	%%%%%%%%%%%

    % Loop over hydraulic conductivity:
    data = struct();
    
    time_factor = 1;
    
    k_sat = 1;
    k_sat_threshold = 1.5;
    
    shock_vector = roundn(exp(linspace(log(1), log(40), 16)), -6);
    data.water_table_elevation_vector = -shock_vector;
    
    num_time_steps = 100;
    data.t_range = zeros(numel(shock_vector), num_time_steps);
    
    data.out_flux = zeros(numel(shock_vector), num_time_steps);
    data.saturation_effective_avg = zeros(numel(shock_vector), 1);

    t_max = zeros(size(shock_vector));
    data.length_vector = zeros(size(shock_vector));
    
    if k_sat < k_sat_threshold
        % Matrix
        t_max = 1.5e+6 * [ ...
            0.000004204123171
            0.000009110296847
            0.000020310761927
            0.000046792350018
            0.000108061512140
            0.000246094813538
            0.000592251355613
            0.001402454295015
            0.003305785123967
            0.007924979825807
            0.019401854479908
            0.046280991735538
            0.110415449258425
            0.272800446048075
            0.646464646464647
            1.561065197428833
            ];
    else
        % Channel
        t_max = 1.5e+5 * [ ...
            0.000001156133872
            0.000002505331633
            0.000005585459530
            0.000012520115275
            0.000028913755951
            0.000067676073723
            0.000158467254610
            0.000385674931129
            0.000909090909091
            0.002179369452097
            0.005191307009489
            0.012727272727273
            0.030364248546067
            0.072992551780431
            0.177777777777778
            0.429292929292929
            ];
    end

    for i = 1:numel(shock_vector)
        max_t = t_max(i);
        data.t_range(i, :) = linspace(0, max_t, num_time_steps);
        dt = data.t_range(i, 2);
        [data.out_flux(i, :), data.saturation_effective_avg(i), t_max(i), vg_par, domain_name] = ...
            run_richards_fdm_implicit(data.water_table_elevation_vector(i), k_sat, data.t_range(i, :), dt, k_sat_threshold, shock_vector(i));
        data.length_vector(i) = shock_vector(i);
        plot(data.t_range(i, :), data.out_flux(i, :));
    end
    
    % t_max .* data.k_sat_vector
    
    data.t_range = data.t_range / time_factor;
    data.van_genuchten_params = vg_par;
    save(sprintf('../Common/data_length_loop_%s_domain.mat', domain_name), '-struct', 'data');
end