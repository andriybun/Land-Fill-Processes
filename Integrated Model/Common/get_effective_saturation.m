function se_3d = get_effective_saturation(spatial_params, hydraulic_params)
    addpath('../Richards/');

    % Calculate hydraulic pressure; we assume water table at the bottom
    % of landfill
    nd = ndims(spatial_params.column_height_array);
    num_domains = size(spatial_params.column_height_array, 2);
    
    ratio = 10;
    fine_dz = spatial_params.dz / ratio;

    % Calculate effective saturation for fine resolution of one column
    hw = fine_dz * cumsum(ones(spatial_params.column_height_array * ratio, 1)) - ...
        spatial_params.zn * spatial_params.dz - fine_dz / 2;
    hw = repmat(hw, [1, num_domains]);
    [~, ~, fine_se] = van_genuchten(hw, hydraulic_params);
    
    % Switch back to coarse resolution
    se = mean(fine_se, 1);
    
    % Repeat a single column over the entire system
    sz = size(spatial_params.column_height_array);
    se_3d = repmat(se, sz) .* (spatial_params.column_height_array > 0);
end