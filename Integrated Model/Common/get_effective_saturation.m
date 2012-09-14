function se_3d = get_effective_saturation(spatial_params, hydraulic_params)
    addpath('../Richards/');

    % Calculate hydraulic pressure; we assume water table at the bottom
    % of landfill
    nd = ndims(spatial_params.is_landfill_array);
    
    ratio = 10;
    fine_dz = spatial_params.dz / ratio;

    % Calculate effective saturation for fine resolution of one column
    hw = fine_dz * cumsum(ones(size(spatial_params.is_landfill_array, 1) * ratio, 1)) - spatial_params.zn - fine_dz / 2;
    [~, ~, fine_se] = van_genuchten(hw, hydraulic_params);
    
    % Switch back to coarse resolution
    se = mean(reshape(fine_se, [], spatial_params.zn))';
    
    % Repeat a single column over the entire system
    sz = size(spatial_params.is_landfill_array);
    se_3d = repmat(se, [1 sz(2:end)]) .* spatial_params.is_landfill_array;
end