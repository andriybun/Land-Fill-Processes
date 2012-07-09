classdef log_normal_params
    properties
        opt_params;         % structure to which a file with best approximations of parameters will be loaded
        se_low;             % minimum value for precomputed effective saturation
        se_hi;              % maximum value for precomputed effective saturation 
    end
    
    methods
        function self = log_normal_params(file_name)
            if nargin == 0
                self.opt_params = load('opt_params_wt.mat');
            else
                self.opt_params = load(file_name);
            end
            if self.opt_params.saturation_effective_avg(1) > self.opt_params.saturation_effective_avg(end)
                self.se_hi = self.opt_params.saturation_effective_avg(1);
                self.se_low = self.opt_params.saturation_effective_avg(end);
            else
                self.se_low = self.opt_params.saturation_effective_avg(1);
                self.se_hi = self.opt_params.saturation_effective_avg(end);
            end
            
        end
        
        function [mu, sigma] = get_params(self, k_sat, se)
            mu = zeros(size(se));
            sigma = zeros(size(se));
%             if k_sat <= 0
%                 error('Hydraulic conductivity must be real positive value');
%             end
%             if (se <= 0) || (se > 1)
%                 error('Effective saturation must be within range 0 to 1');
%             end
            se = min(se, self.se_hi);
            se = max(se, self.se_low);

            mu = - log(k_sat) + interp1(self.opt_params.saturation_effective_avg, self.opt_params.z, se);
            sigma = interp1(self.opt_params.saturation_effective_avg, self.opt_params.v, se);
        end
    end
end