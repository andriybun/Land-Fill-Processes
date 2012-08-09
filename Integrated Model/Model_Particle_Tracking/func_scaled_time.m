function y = func_scaled_time(t, u, f_handle)
    y = nan(size(t));
    udt = 0;
    u0 = u(1);

    % step t = 0
    y(1) = f_handle(t(1), u0);

    % loop over t
    for idx = 2:numel(t)
        ratio = u(idx) / u0;
        udt = udt + (t(idx) - t(idx-1)) * ratio;
        y(idx) = ratio * f_handle(udt, u0);
    end
end