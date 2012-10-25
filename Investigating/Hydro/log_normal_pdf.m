function res = log_normal_pdf(t_in, mu_in, sigma_in, dt)
    t_in = repmat(t_in, [numel(mu_in), 1]);
    mu_in = repmat(mu_in, [1, size(t_in, 2)]);
    sigma_in = repmat(sigma_in, [1, size(t_in, 2)]);
    res = exp(-(log(t_in) - mu_in).^2 ./ (2 .* sigma_in .* sigma_in)) ./ (sqrt(2 * pi) .* sigma_in .* t_in);
    sigma_is_inf = isinf(sigma_in(:, 1));
    res(~sigma_is_inf, 1) = 0;
    res(sigma_is_inf, 1) = 1 / dt;
end
