function properties_array = generate_biogeochemical_properties(geometry_params, hydraulic_params)
    if isfield(geometry_params, 'num_domains')
        is_landfill_array = repmat(geometry_params.is_landfill_array, [1, geometry_params.num_domains]);
    else
        is_landfill_array = geometry_params.is_landfill_array;
    end
    
    is_data_idx = (is_landfill_array > 0);
    sz = size(is_landfill_array);
    
    properties_array_soa = struct();
    
    num_dims = isfield(geometry_params, 'dz') + isfield(geometry_params, 'dx') + isfield(geometry_params, 'dy');
    
    if num_dims == 1
        properties_array_soa.effective_saturation = get_effective_saturation_1d(geometry_params, hydraulic_params);
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