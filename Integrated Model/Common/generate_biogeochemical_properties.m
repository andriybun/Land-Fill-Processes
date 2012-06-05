function properties_array = generate_biogeochemical_properties(spatial, column_height_array)
    sz = size(column_height_array);
    properties_array_soa = struct();
    
    properties_array_soa.effective_saturation = 0.3 * ones(sz);
    
    %% STUB:
    properties_array_soa.compound_1_fraction = zeros(sz);
    properties_array_soa.compound_1_fraction(column_height_array > 0) = 0.5 * ones(spatial.num_columns, 1);
    properties_array_soa.compound_2_fraction = zeros(sz);
    properties_array_soa.compound_2_fraction(column_height_array > 0) = 0.3 * ones(spatial.num_columns, 1);
    %% END STUB
    
    % Convert SoA to AoS
    properties_array = soa_to_aos(properties_array_soa);
end