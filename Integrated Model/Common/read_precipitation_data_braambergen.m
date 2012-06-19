function [precipitation_intensity, time_params, start_date] = read_precipitation_data_braambergen(file_name)
    if nargin < 1
        file_name = '../Data/precip_braambergen.mat';
    end
    load(file_name);
end