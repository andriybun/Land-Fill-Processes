function [f, q_in, dhw_above] = net_flux(hwn, t, dt, theta_prev, k_n, in_flx, hw_above, spatial, vg_par, h_out, accumulate_water, EPSILON)
    nn = spatial.nn;
    % Computing internode hydraulic head:
    h_in = zeros(nn+1, 1);
    k_in = zeros(nn+1, 1);
    h_in(2:end-1) = (hwn(1:end-1) + hwn(2:end)) ./ 2;
    % Computing internode conductivity
    [k_in(2:end-1), ~] = van_genuchten(h_in(2:end-1), vg_par);
    % Internode conductivity at system boundaries = maximum
    k_in(1) = k_n(1);
    k_in(nn+1) = k_n(nn);
    % Append hydraulic head below the system to the vector
    hwn = vertcat(hwn, h_out);
    % Compute moisture content
    [~, theta_curr] = van_genuchten(hwn, vg_par);
    % Compute internode flux
%         in_flux_limit = -spatial.dzin(1) / dt * (vg_par.theta_s - theta_curr(1));
    in_flux_limit = k_in(1);
    dhw_above = 0;
    q_in = flux(hwn, k_in, spatial, in_flx, in_flux_limit);
    f = (theta_prev(1:nn) - theta_curr(1:nn)) ./ dt - (q_in(2:nn+1) - q_in(1:nn)) ./ spatial.dzin(1:nn);

    return;
    
    function qx = flux(h_nx, k_inx, spatial, in_flx, in_flux_limit)
        nn = spatial.nn;
        dzx = vertcat(spatial.dzn(1), spatial.dzn, spatial.dzn(nn-1));
        qx = zeros(nn+1, 1);
        qx(1) = in_flx - hw_above;
        dhw_above = - hw_above;
        hw_above = 0;
        qx(2:nn+1) = - k_inx(2:nn+1) .* ((h_nx(2:nn+1) - h_nx(1:nn)) ./ dzx(2:nn+1) + 1);
        if accumulate_water
            flux_potential = real_lt(in_flux_limit + qx(1), 0, EPSILON);
            if flux_potential || ~real_eq(hw_above, 0, EPSILON)
                qx(1) = -min(in_flux_limit, -qx(1));
                in_flux_limit_remaining = in_flux_limit + qx(1);
                if (real_lt(0, in_flux_limit_remaining, EPSILON))
                    dhw_above = - min(in_flux_limit_remaining, hw_above);
                    qx(1) = qx(1) + dhw_above;
                end
                dhw_above = qx(1) - in_flx;
            end
        else
            qx(1) = max(-k_inx(1), qx(1));
        end
%         qx(1) = max(-k_inx(1), qx(1));
%         qx(1) = min(k_inx(1), qx(1));
        qx(nn+1) = max(-k_inx(nn+1), qx(nn+1));
        qx(nn+1) = min(k_inx(nn+1), qx(nn+1));
    end

end