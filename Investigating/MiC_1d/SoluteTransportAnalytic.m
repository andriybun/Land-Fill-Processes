function cOut = soluteTransportAnalytic(x, t, u, d)
    nx = numel(x);
    nt = numel(t);

    checkInputs(x, u);
    checkInputs(x, d);
    
    u = -u;
    
    if ~(nx == 1 || nt == 1)
        x = repmat(x, [1, nt]);
        t = repmat(t, [nx, 1]);
        if numel(u) ~= 1
            u = repmat(u, [1, nt]);
        end
        if numel(d) ~= 1
            d = repmat(d, [1, nt]);
        end
    end
    
    x = -x;
    l = x(end);
    %% IC: C = 1, 0 <= x <= l
    % cOut = -0.5 .* erf(0.5 .* (-x + u .* t) ./ sqrt(d .* t)) ./ l + ...
    %     0.5 .* erf(0.5 .* (-x + l + u .* t) ./ sqrt(d .* t)) ./ l;
    %% IC: C = 1, 0 <= x < inf
    cOut = 0.5 - 0.5 * erf(0.5 * (-x + u .* t) ./ (sqrt(t .* d)));

    %%
    return
    
    function checkInputs(spatial, var)
        if (numel(var) ~= 1) && (~isequal(size(var), size(spatial)))
            error('Wrong dimensions of input arguments');
        end
    
    end
end