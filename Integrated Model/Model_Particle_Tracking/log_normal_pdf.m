function res = log_normal_pdf(t, mu, sigma)
    if iscolumn(mu)
        % 2D case
        t = repmat(t, [numel(mu), 1]);
        mu = repmat(mu, [1, size(t, 2)]);
        sigma = repmat(sigma, [1, size(t, 2)]);
%     else
%         % 3D case
%         t = repmat(t', [1, size(mu)]);
%         t = permute(t, [2, 3, 1]);
%         mu = repmat(mu, [1, 1, size(t, 3)]);
%         sigma = repmat(sigma, [1, 1, size(t, 3)]);
    end
    res = exp(-(log(t) - mu).^2 ./ (2 .* sigma .* sigma)) ./ (sqrt(2 * pi) .* sigma .* t);
end