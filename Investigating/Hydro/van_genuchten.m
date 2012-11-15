function [k, theta, se] = van_genuchten(hw, vg_par)
    alpha   = vg_par.alpha;
    theta_s = vg_par.theta_s;
    theta_r = vg_par.theta_r;
    n       = vg_par.n;
    m       = vg_par.m;
    k_sat   = vg_par.k_sat;
    
    % Calculate the capillary pressure
    hc = 0 - hw;
    
    % Compute the effective saturation
    se = (1 + (alpha .* hc).^n).^(-m) .* (hc > 0) + 1 .* (hc <= 0);
    
    % Compute the effective saturation
    theta = se .* (theta_s - theta_r) + theta_r;
 
    % Compute the hydraulic conductivity
    k = k_sat .* se.^(1/2) .* (1 - (1 - se.^(1./m)).^m).^2;
 
%     % Compute the specific moisture storage
%     c = -alpha .* n .* (hc) .* (1 / n - 1) .* (alpha .* (hc)).^(n - 1) .*...
%        (theta_r - theta_s) .* ((alpha .* (hc)).^n + 1).^(1/n - 2) .* (hc > 0) + 0 * (hc <= 0);
end