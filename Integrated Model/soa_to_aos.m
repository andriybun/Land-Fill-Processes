function aos = soa_to_aos(soa)
    % Function to convert structure of arrays to an array of structures
    
    % Make sure all arrays are the same sizes, and generate a
    % zero-structure with the same fields
    [zero_struct, sz] = preprocess(soa);
    % Initialize array of zero-structures
    aos = repmat(zero_struct, sz);
    % Fill in array of structures with original data
    aos = fill_aos(soa, aos, sz);
    
    return
    
    function [zero_struct, sz] = preprocess(soa)
        names = fieldnames(soa);
        zero_struct = struct();
        if ~isstruct(soa.(names{1}))
            sz = size(soa.(names{1}));
            zero_struct.(names{1}) = 0;
        else
            [zero_struct.(names{1}), sz] = preprocess(soa.(names{1}));
        end
        for i = 2:numel(names)
            name = names{i};
            if ~isstruct(soa.(name))
                zero_struct.(name) = 0;
            else
                [zero_struct.(name), ~] = preprocess(soa.(name));
            end
            if (~isempty(find(((size(soa.(name)) ~= sz) ~= 0))))
                error(sprintf('Error! Inconsistent array sizes: size(%s) ~= size(%s)', names{1}, name));
            end
        end
    end

    function aos = fill_aos(soa, aos, sz)
        names = fieldnames(soa);
        num_elements = prod(sz);
        if ~isstruct(soa.(names{1}))
            for idx = 1:num_elements
                aos(idx).(names{1}) = soa.(names{1})(idx);
            end
        else
            internal_aos = soa_to_aos(soa.(names{1}));
            for idx = 1:num_elements
                aos(idx).(names{1}) = internal_aos(idx);
            end
        end
        for i = 2:numel(names)
            name = names{i};
            if ~isstruct(soa.(name))
                for idx = 1:num_elements
                    aos(idx).(name) = soa.(name)(idx);
                end
            else
                internal_aos = soa_to_aos(soa.(name));
                for idx = 1:num_elements
                    aos(idx).(name) = internal_aos(idx);
                end
            end
        end
    end
end