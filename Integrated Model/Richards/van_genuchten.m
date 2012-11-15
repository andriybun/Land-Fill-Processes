function [k, theta, se] = van_genuchten(hw, vg_par)
 
    alpha   = reproduce(hw, vg_par.alpha);
    theta_s = reproduce(hw, vg_par.theta_s);
    theta_r = reproduce(hw, vg_par.theta_r);
    n       = reproduce(hw, vg_par.n);
    m       = reproduce(hw, vg_par.m);
z    k_sat   = reproduce(hw, vg_par.k_sat);
    
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
 
    return
    
    function res_ar = reproduce(sample, ar)
        sz = size(sample);
        if numel(ar) == 1 || isequal(sz, size(ar))
            res_ar = ar;
        elseif sz(1) == size(ar, 1)
            res_ar = repmat(ar, [1, sz(2)]);
        elseif sz(2) == size(ar, 2)
            res_ar = repmat(ar, [sz(1), 1]);
        else
            err = MException('InputChk:IncorrectSize', 'Error! Array sizes are inconsistent.');
            throw(err);
        end
    end
endz