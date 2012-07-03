function properties_array = generate_biogeochemical_properties_3d(spatial_params)
    is_landfill_array = spatial_params.is_landfill_array;
    sz = size(is_landfill_array);
    properties_array_soa = struct();
    
    properties_array_soa.effective_saturation = 0.5 * ones(sz);
    
    %% STUB:
    is_data_idx = (is_landfill_array > 0);
    properties_array_soa.solutes.solute_1_fraction = zeros(sz);
    properties_array_soa.solutes.solute_1_fraction(is_data_idx) = 0.5;
    properties_array_soa.solutes.solute_2_fraction = zeros(sz);
    properties_array_soa.solutes.solute_2_fraction(is_data_idx) = 0.3;
    %% END STUB
    
%     % Convert SoA to AoS
%     properties_array = soa_to_aos(properties_array_soa);
    properties_array = properties_array_soa;
end