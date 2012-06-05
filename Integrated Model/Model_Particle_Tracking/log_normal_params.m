classdef log_normal_params
    properties
        opt_params;         % structure to which a file with best approximations of parameters will be loaded
    end
    
    methods
        function self = log_normal_params()
            self.opt_params = load('opt_params_wt.mat');
        end
        
        function [mu, sigma] = get_params(self, k_sat, se)
%             if k_sat <= 0
%                 error('Hydraulic conductivity must be real positive value');
%             end
%             if (se <= 0) || (se > 1)
%                 error('Effective saturation must be within range 0 to 1');
%             end
            mu = - log(k_sat) + interp1(self.opt_params.saturation_effective_avg, self.opt_params.z, se);
            sigma = interp1(self.opt_params.saturation_effective_avg, self.opt_params.v, se);
        end
    end
end