function [mu, sigma] = compute_lognormal_parameters(mean, confidence_interval)
    
    sigma = fzero(@confidence_difference, 100);
    mu = log(mean) - sigma^2 / 2;

    return
    
    function diff = confidence_difference(sigma)
        diff = confidence_interval - mean * exp(-1 / 2 * sigma + 1.96 * sigma) + mean * exp(-1 / 2 * sigma - 1.96 * sigma);
    end
end