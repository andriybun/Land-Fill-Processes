function [precipitation_intensity, time_params, start_date] = read_precipitation_data_csv_braambergen(file_name)
    if nargin < 1
        file_name = '../Data/MeteoAll.csv';
    end

    [raw, delimiter, nheaderlines] = importdata(file_name);

    header = {raw.textdata{nheaderlines, :}};

    [nr, nc] = size(raw.data);

    % Parse time data
    time = nan(nr, 6);
    time_num = nan(nr, 1);
    for idx = 1:nr,
        time(idx, :) = datevec(raw.textdata{idx + nheaderlines, 1}, 'dd-mm-yy HH:MM:SS');
        time_num(idx) = datenum(time(idx, 1:3));
    end

    % Remove data for non exact hours
    drop_idx = find(time(:,5) ~= 0);
    time(drop_idx, :) = [];
    time(:, 5:6) = [];
    time_num(drop_idx) = [];
    data = raw.data;
    data(drop_idx, :) = [];
    
    % Start date
    start_date = struct();
    start_date.year = time(end, 1);
    start_date.month = time(end, 2);
    start_date.day = time(end, 3);
    
%     % Validate dataset
%     val = time(1:end-1, :) - time(2:end, :);
%     for idx = size(val, 1):-1:1
%         if val(idx, 4) ~= 1 || sum(abs(val(idx, 1:3))) ~= 0
%             if val(idx, 4) == -23
%                 if (val(idx, 3) ~= 1) && (val(idx, 2) ~= 1) && (val(idx, 1) ~= 1)
%                     disp(time(idx:idx+1, :));
%                 end
%             else
%                 disp(time(idx:idx+1, :));
%             end
%         end
%     end
    
    % Time params
    time_params.max_days = time_num(1)-time_num(end)+1;                                             % number of simulation days
    time_params.time_discretization = 3600;                                                         % 1 hr in seconds
    time_params.intervals_per_day = 24 * 3600 / time_params.time_discretization;
    num_intervals = time_params.max_days * time_params.intervals_per_day;                           % in {time step}
    time_params.num_intervals = num_intervals;
    time_params.days_elapsed = (0 : 1: (num_intervals-1)) / time_params.intervals_per_day;
    
    % Prepare consistent time array
    time_processed = nan(num_intervals, 4);
    time_processed(:, 4) = mod(0:num_intervals-1, 24);
    tmp_datenum = reshape(repmat((time_num(end):time_num(1)),[24, 1]), [], 1);
    tmp_datevec = datevec(tmp_datenum);
    time_processed(:, 1:3) = tmp_datevec(:, 1:3);

    % Parse precipitation data
    precipitation_cumulative = zeros(num_intervals, 1);
    
    raw_idx = numel(time_num);
    prev_precipitation = data(raw_idx, 2);
    for idx = 1:num_intervals
        if (raw_idx < 1) || ~date_eq(time_processed(idx, :), time(raw_idx, :))
            while date_lt(time(raw_idx, :), time_processed(idx, :))
                raw_idx = raw_idx - 1;
            end
            precipitation_cumulative(idx) = prev_precipitation;
            continue
        end
        precipitation_cumulative(idx) = data(raw_idx, 2);
        prev_precipitation = precipitation_cumulative(idx);
        raw_idx = raw_idx - 1;
    end
    
    precipitation_cumulative = precipitation_cumulative * 1e-4;
    precipitation_intensity = precipitation_cumulative - vertcat(0, precipitation_cumulative(1:end-1));
    
    save('../Data/precip_braambergen.mat');
    
    return
    
    function res = date_eq(d1, d2)
        res = (sum(abs(d1 - d2)) == 0);
    end

    function res = date_lt(d1, d2)
        res = (date_time_key(d1) < date_time_key(d2));
    end

    function res = date_time_key(d)
        res = 24 * 366 * d(:, 1) + 24 * 31 * d(:, 2) + 24 * d(:, 3) + d(:, 4);
    end
end