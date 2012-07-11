function res = log_normal_cdf(t, mu, sigma)
    t = repmat(t, [numel(mu), 1]);
    mu = repmat(mu, [1, size(t, 2)]);
    sigma = repmat(sigma, [1, size(t, 2)]);
    res = 1 / 2 + 1 / 2 * erf((log(t) - mu) ./ sqrt(2 .* sigma .^ 2));
    sigma_is_inf = isinf(sigma(:, 1));
    res(sigma_is_inf, 1) = 0;
    res(sigma_is_inf, 2:end) = 1;
end