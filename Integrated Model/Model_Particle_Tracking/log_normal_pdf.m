function res = log_normal_pdf(t, mu, sigma)
    t = repmat(t, [numel(mu), 1]);
    mu = repmat(mu, [1, size(t, 2)]);
    sigma = repmat(sigma, [1, size(t, 2)]);
    res = exp(-(log(t) - mu).^2 ./ (2 .* sigma .* sigma)) ./ (sqrt(2 * pi) .* sigma .* t);
end