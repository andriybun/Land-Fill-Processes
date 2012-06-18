function precipitation_intensity = generate_precipitation_data(start_date, time_params)
    
    rainy_time_fraction = 0.06;
    
    num_intervals = time_params.num_intervals;
    intervals_per_day = time_params.intervals_per_day;
    max_days = time_params.max_days;
    
    % Initialize some arrays
    interval_months = zeros(1, num_intervals);
    is_precipitation = (rand(1, num_intervals) <= rainy_time_fraction); % Occurrence of a rainfall event during interval
    
%     %% STUB: single input pulse
%     is_precipitation = zeros(1, num_intervals);
%     is_precipitation(1:1) = 1;
%     %% END STUB
    
    precipitation_intensity = zeros(1, num_intervals);
    
    % Number of days and then intervals in the first month
    current_month = start_date.month;
    first_month_days_remaining = month_length(current_month) - start_date.day + 1;
    interval_months(1:min(num_intervals, first_month_days_remaining * intervals_per_day)) = current_month;
    current_day = first_month_days_remaining;                   % Last day of the last accounted month
    days_remaining = max_days - first_month_days_remaining;     % Number of days for calculation remaining
    
    % Remaining months
    while (days_remaining > 0)
        current_month = (current_month < 12) * (current_month + 1) + (current_month == 12) * 1;
        days_current_month = month_length(current_month);
        last_interval = min(num_intervals, (current_day + days_current_month) * intervals_per_day);
        interval_months((current_day * intervals_per_day + 1) : last_interval) = current_month;
        current_day = current_day + days_current_month;
        days_remaining = days_remaining - days_current_month;
    end
    
    precipitation_intensity = precipitation_statistics(interval_months) ./ month_length(interval_months) / ...
        intervals_per_day / rainy_time_fraction;
    precipitation_intensity = precipitation_intensity .* is_precipitation * 1e-3;
    
    return
    
    %% Precipitation statistics by months
    function precip_res = precipitation_statistics(months)
        % Precipitation statistics by months:
        % Precipitation in mm
        precipitation = [ 68 53 44 49 52 58 77 87 72 72 70 64 ];
        % Days with precipitation
        precip_days = [ 22 19 16 16 14 14 17 18 19 20 21 21 ];
        precip_res = precipitation(months);
    end

    function len = month_length(month)
        eom_day = [31 28 31 30 31 30 31 31 30 31 30 31];
        len = eom_day(month);
    end

%     %% Function stub:
%     function [start_date, max_days, time_discretization] = stub()
%         % Time
%         start_date = struct();
%         start_date.month = 1;
%         start_date.day = 1;
%         max_days = 365;                                         % number of simulation days
%         time_discretization = 3600;                             % in seconds
%         max_time = max_days * 24 * 3600 / time_discretization;  % in {time step}
%     end
end