function res = log_normal_pdf(t, mu, sigma, dt)
    t = repmat(t, [numel(mu), 1]);
    mu = repmat(mu, [1, size(t, 2)]);
    sigma = repmat(sigma, [1, size(t, 2)]);
    res = exp(-(log(t) - mu).^2 ./ (2 .* sigma .* sigma)) ./ (sqrt(2 * pi) .* sigma .* t);
    sigma_is_inf = isinf(sigma(:, 1));
    res(~sigma_is_inf, 1) = 0;
    res(sigma_is_inf, 1) = 1 / dt;
end