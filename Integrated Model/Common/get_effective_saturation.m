function se = get_effective_saturation(spatial_params, hydraulic_params)
    addpath('../Richards/');

    % Calculate hydraulic pressure; we assume water table at the bottom
    % of landfill
    nd = ndims(spatial_params.is_landfill_array);
    if nd == 2
        hw = cumsum(ones(size(spatial_params.is_landfill_array)), 1) - spatial_params.zn - spatial_params.dz / 2;
    else
        hw = cumsum(ones(size(spatial_params.is_landfill_array)), nd) - spatial_params.zn - spatial_params.dz / 2;
    end
    hw = hw .* spatial_params.is_landfill_array;
    [~, ~, se] = van_genuchten(hw, hydraulic_params);
end