function saturation_effective_avg = calc_avg_se(water_table_elevation, zn, vg_par)
    hw = water_table_elevation - zn;
    
    [~, theta(:, 1)] = van_genuchten(hw, vg_par);

    % Total relative moisture content
    saturation_effective_avg = mean((theta(:, 1) - vg_par.theta_r) / (vg_par.theta_s - vg_par.theta_r));
end