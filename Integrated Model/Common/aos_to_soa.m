function soa = aos_to_soa(aos)
    % Converting array of structures to structure of arrays
    % Non recursive yet!!!
    names = fieldnames(aos(1));
    num_arrays = numel(aos);
    for i = 1:numel(names)
        name = names{i};
        soa.(name) = zeros(size(aos));
        for idx = 1:num_arrays
            soa.(name)(idx) = aos(idx).(name);
        end
    end
end