function properties_array = generate_biogeochemical_properties_3d(spatial_params, hydraulic_params)
    addpath('../Richards/');

    is_landfill_array = spatial_params.is_landfill_array;
    is_data_idx = (is_landfill_array > 0);
    sz = size(is_landfill_array);
    
    properties_array_soa = struct();
    if nargin < 2
        properties_array_soa.effective_saturation = 0.3 * ones(sz);
    else
        % Calculate hydraulic pressure; we assume water table at the bottom
        % of landfill
        hw = cumsum(ones(size(spatial_params.is_landfill_array)), 3) - spatial_params.zn - spatial_params.dz / 2;
        hw = hw .* spatial_params.is_landfill_array;
        [~, ~, properties_array_soa.effective_saturation] = van_genuchten(hw, hydraulic_params);
    end
    
    %% STUB:
    properties_array_soa.solutes.solute_1_fraction = zeros(sz);
    properties_array_soa.solutes.solute_1_fraction(is_data_idx) = 0.5;
    properties_array_soa.solutes.solute_2_fraction = zeros(sz);
    properties_array_soa.solutes.solute_2_fraction(is_data_idx) = 0.3;
    %% END STUB
    
%     % Convert SoA to AoS
%     properties_array = soa_to_aos(properties_array_soa);
    properties_array = properties_array_soa;
end