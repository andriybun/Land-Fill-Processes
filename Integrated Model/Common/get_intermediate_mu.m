function mu_steps = get_intermediate_mu(geometry_params, log_normal_params, mean_steps_rel)
% Calculating log-normal parameters for intermediate points

    sz = size(geometry_params.is_landfill_array);
    if numel(sz) == 2
        sz = [sz, 1];
    end
    
    % Average flux time at the bottom (baseline)
    mean_orig = exp(log_normal_params.mu + log_normal_params.sigma .* log_normal_params.sigma / 2);
    % Values of average flux time for intermediate points
    mean_steps = permute(repmat(mean_orig, [1, 1, 1, numel(mean_steps_rel)]), [4, 1, 2, 3]) .* repmat(mean_steps_rel', [1, size(mean_orig)]);

    % Mu for corresponding distributions of out-flux at intermediate points
    mu_steps = log(mean_steps(2:end, :, :, :)) - ...
        permute(repmat(log_normal_params.sigma .* log_normal_params.sigma / 2, [1, 1, 1, geometry_params.zn]), [4, 1, 2, 3]);
    
    for x_idx = 1:sz(2)
        for y_idx = 1:sz(3)
            z_idx = 1;
            while (z_idx <= geometry_params.zn)
                z_idx_start = z_idx;
                shift_flag = false;
                while (z_idx <= geometry_params.zn) && (geometry_params.is_landfill_array(z_idx, x_idx, y_idx) == 0)
                    z_idx = z_idx + 1;
                    shift_flag = true;
                end
                if shift_flag
                    num_shift = geometry_params.zn - z_idx + 1;
                    mu_steps(z_idx:end, :, :, :) = mu_steps(z_idx_start:z_idx_start + num_shift - 1, :, :, :);
                    if z_idx_start == 1
                        mu_steps(z_idx_start:z_idx-1, :, :, :) = 0;
                    else
                        mu_steps(z_idx_start:z_idx-1, :, :, :) = mu_steps(z_idx_start-1, :, :, :);
                    end
                end
                z_idx = z_idx + 1;
            end
        end
    end
end