function se = get_effective_saturation_1d(geometry_params, hydraulic_params)
    addpath('../Richards/');

    % Calculate hydraulic pressure; we assume water table at the bottom
    % of landfill
 
    if isfield(geometry_params, 'num_domains')
        num_domains = geometry_params.num_domains;
    else
        num_domains = 1;
    end
    
    ratio = 10;
    fine_dz = geometry_params.dz / ratio;

    % Calculate effective saturation for fine resolution of one column
    hw = fine_dz * cumsum(ones(geometry_params.column_height_array * ratio, 1)) - ...
        geometry_params.zn * geometry_params.dz - fine_dz / 2;
    hw = repmat(hw, [1, num_domains]);
    [~, ~, fine_se] = van_genuchten(hw, hydraulic_params);
    
    % Switch back to coarse resolution
    se = mean(fine_se, 1);
end