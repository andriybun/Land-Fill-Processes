function [precipitation_intensity, time_params, start_date] = read_precipitation_data_csv(file_name)
    if nargin < 1
        file_name = '../Data/precipitation_daily_KNMI_20110908.txt';
    end

    % Read daily precipitation data
    idx = 1;
    fid = fopen(file_name);
    while (~feof(fid))
        raw_line = fgets(fid);
        if raw_line(1) ~= '#'
            raw_data = textscan(raw_line, '%*s %08s %f %f %f', 'delimiter', ',');
            if idx == 1
                start_date = struct();
                start_date.year = str2double(raw_data{1}{1}(1:4));
                start_date.month = str2double(raw_data{1}{1}(5:6));
                start_date.day = str2double(raw_data{1}{1}(7:8));
            end
            [out.date(idx), out.precipitation_duration(idx), out.precipitation(idx), out.evaporation(idx)] = parse_data(raw_data);
            idx = idx + 1;
        end
    end
    fclose(fid);
    
    % Adjust precipitation to evaporation data
    out.precipitation = out.precipitation * ((sum(out.precipitation) - sum(out.evaporation)) / sum(out.precipitation));
    
    % Adjust time_params
    time_params.max_days = 60; % numel(out.precipitation);  % number of simulation days
    time_params.time_discretization =     3600;                                                     % 0.1 hrs in seconds
    time_params.intervals_per_day = 24 * 3600 / time_params.time_discretization;
    num_intervals = time_params.max_days * time_params.intervals_per_day;                           % in {time step}
    time_params.num_intervals = num_intervals;
    time_params.days_elapsed = (0 : 1: (num_intervals-1)) / time_params.intervals_per_day;
    
    % Generate precipitation intensity vector
    precipitation_intensity = zeros(1, num_intervals);
    
    for i = 1:time_params.max_days
        start_of_day_idx = time_params.intervals_per_day * (i - 1) + 1;
        precipitation_intensity(start_of_day_idx:start_of_day_idx + out.precipitation_duration(i) - 1) = ...
            out.precipitation(i) / out.precipitation_duration(i);
    end
    
    return
    
    function [date, precipitation_duration, precipitation, evaporation] = parse_data(raw_data)
        date = datenum(str2double(raw_data{1}{1}(1:4)), str2double(raw_data{1}{1}(5:6)), str2double(raw_data{1}{1}(7:8)));
        precipitation_duration = ceil(raw_data{2} / 10);   % hrs
        precipitation = raw_data{3} * 1e-4;                   % m
        evaporation = raw_data{4} * 1e-4;                     % m
    end
end