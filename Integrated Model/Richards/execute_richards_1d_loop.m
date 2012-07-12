function execute_richards_1d_loop()

    %%%%%%%%%%%
    %%%
    %%% first time step is void
    %%%
	%%%%%%%%%%%

    %% Loop over water table:
    data = struct();
    
    time_factor = 1;
    k_sat = 0.1 / time_factor; % m/day
    
    data.water_table_elevation_vector = [-0.1:-0.4:-19.9];
    num_time_steps = 100;
    data.t_range = zeros(numel(data.water_table_elevation_vector), num_time_steps);
    
    data.out_flux = zeros(numel(data.water_table_elevation_vector), num_time_steps);
    data.saturation_effective_avg = zeros(numel(data.water_table_elevation_vector), 1);

    if k_sat < 10
        % Matrix domain:
        t_max = 1e+4 * [ ...
            0.000049469287302
            0.000620614695245
            0.002963941105689
            0.008826520110156
            0.018554886564644
            0.031732299200368
            0.048332992763524
            0.067542066929885
            0.089611365343897
            0.114491418718260
            0.140665918443696
            0.170642057944838
            0.201998590437618
            0.236523156076820
            0.271518929169819
            0.309311965599358
            0.345230978004999
            0.389044464175944
            0.428849938695873
            0.472206833338814
            0.520608033756043
            0.568649457538346
            0.610543761263385
            0.668321077036328
            0.721719893464892
            0.774510899748584
            0.826173075609590
            0.884390931417673
            0.932490857941415
            1.003439520844981
            1.050287525137648
            1.111970973104500
            1.176468387579499
            1.239619594089299
            1.317588233580109
            1.381326412130319
            1.437336257517872
            1.526894287224866
            1.583050740029620
            1.643113216762164
            1.722888660551484
            1.797447440990308
            1.887119643563791
            1.968657771126907
            2.014889693326959
            2.092438839091160
            2.180875698203242
            2.261930094837683
            2.341232173138008
            2.418420413801218
            ];
    else
        % Channel domain:
        t_max = [ ...
            0.000004784066879
            0.000059546364345
            0.000279257180480
            0.000857467646571
            0.001753411405488
            0.003017626525432
            0.004588276397275
            0.006406272705630
            0.008588417934403
            0.010878662716911
            0.013486997052544
            0.016184396463053
            0.019238056173063
            0.022520740361323
            0.025574400071333
            0.029900417993848
            0.033590256810110
            0.037407331447622
            0.041860585191387
            0.045423188186398
            0.049876441930163
            0.053947988210176
            0.059037421060193
            0.064126853910209
            0.068707343475224
            0.073796776325241
            0.078886209175257
            0.083975642025274
            0.090082961445294
            0.095274182952311
            0.100898006251580
            0.107183455821350
            0.113290775241370
            0.119092728690389
            0.125963463037911
            0.131918099472431
            0.139246882776455
            0.145354202196475
            0.151156155645493
            0.159401036862520
            0.165355673297040
            0.171004943760558
            0.179555190948586
            0.185967876339607
            0.194823489498636
            0.200167393991153
            0.207190811324176
            0.216199107468706
            0.225207403613235
            0.234215699757765
            ];
    end;
    for i = 1:numel(data.water_table_elevation_vector)
%         if (data.water_table_elevation_vector(i) >= -1)
%             max_t = 100;
%         elseif (data.water_table_elevation_vector(i) >= -3)
%             max_t = 900;
%         else
%             max_t = 100000;
%         end
        max_t = t_max(i);
        
        data.t_range(i, :) = linspace(0, max_t, num_time_steps);
        dt = data.t_range(i, 2);
        [data.out_flux(i, :), data.saturation_effective_avg(i), t_max(i), vg_par, domain_name] = ...
            run_richards_fdm_implicit(data.water_table_elevation_vector(i), k_sat, data.t_range(i, :), dt);
        plot(data.t_range(i, :), data.out_flux(i, :));
    end
    
    data.t_range = data.t_range / time_factor;
    data.van_genuchten_params = vg_par;
    save(sprintf('../Common/data_wt_loop_%s_domain.mat', domain_name), '-struct', 'data');
    
%     disp([data.water_table_elevation_vector t_max]);
    
%     %% Loop over hydraulic conductivity:
%     data = struct();
%     
%     water_table_elevation = -2;
% 
%     data.k_sat_vector =  1e+3 * [0.001, 0.003, 0.007, 0.01, 0.03, 0.07, 0.1, 0.3, 0.7, 1, 5, 10];
%     
%     num_time_steps = 100;
%     data.t_range = zeros(numel(data.k_sat_vector), num_time_steps);
%     data.out_flux = zeros(numel(data.k_sat_vector), num_time_steps);
%     
%     for i = 1:numel(data.k_sat_vector)
%         max_t = 35 / data.k_sat_vector(i);
%         data.t_range(i, :) = linspace(0, max_t, num_time_steps);
%         dt = data.t_range(i, 2);
%         [data.out_flux(i, :), ~] = run_richards_fdm_implicit(water_table_elevation, data.k_sat_vector(i), data.t_range(i, :), dt);
%     end
%     
%     save('data_hydro_conductivity_loop.mat', '-struct', 'data');

end