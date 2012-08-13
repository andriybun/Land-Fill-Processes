% Solute residence time generator class

classdef solute_residence_generator
    
    properties (Access = public)
        l;
        d;
        u;
    end
    
    properties (Access = private)
        EPSILON;
        cdf_handle;
        sz;
        batch_size = 50;
        num_batches;
        l_is_vector;
        d_is_vector;
        u_is_vector;
    end
    
    methods (Access = public)
        function self = solute_residence_generator(l, d, u, EPSILON, ic_type, rng_seed)
            self.l = l;
            self.d = d;
            self.u = u;
            self.EPSILON = EPSILON;
            
            % If passed, set the random seed
            if nargin > 5
                rng(rng_seed);
            end
            
            % Validate inputs' dimensions
            try
                l + u + d;
            catch exception
                error('Dimensions error: input parameters'' dimensions do not agree')
            end
            
            % Set internal size of domain
            self.sz = [1, 1];
            if numel(l) > 1
                self.sz = size(l);
                self.l_is_vector = true;
            elseif numel(u) > 1
                self.sz = size(u);
                self.u_is_vector = true;
            elseif numel(d) > 1
                self.sz = size(d);
                self.d_is_vector = true;
            end
            
            self.num_batches = ceil(prod(self.sz) / self.batch_size);
            
            switch ic_type
                case 1
                    self.cdf_handle = @self.cdf_ic;
                case 2;
                    self.cdf_handle = @self.cdf_if;
                otherwise
                    error('Error! Wrong IC type.');
            end
        end
        
        function t_gen = generate(self)
            rn = unifrnd(0, 1, self.sz);
            num_cells = prod(self.sz);
            ini_t = ones(self.sz);
            options = optimset('Display', 'off'); % , 'TolX', self.EPSILON
            t_gen = nan(self.sz);
            first_cell_idx = 1;
            for idx = 1:self.num_batches
                last_cell_idx = min(num_cells, first_cell_idx + self.batch_size - 1);
                t_gen(first_cell_idx:last_cell_idx) = lsqnonlin(@f, ini_t(first_cell_idx:last_cell_idx), ...
                    zeros(last_cell_idx - first_cell_idx + 1, 1), inf(last_cell_idx - first_cell_idx + 1, 1) , options);
                first_cell_idx = last_cell_idx + 1;
            end
            return
            
            function y = f(t)
                l = self.l;
                d = self.d;
                u = self.u;
                if self.l_is_vector
                    l = l(first_cell_idx:last_cell_idx);
                end
                if self.d_is_vector
                    d = d(first_cell_idx:last_cell_idx);
                end
                if self.u_is_vector
                    u = u(first_cell_idx:last_cell_idx);
                end
                y = self.cdf_handle(t, l, d, u) - rn(first_cell_idx:last_cell_idx); 
            end
        end
        
        function self = set_l(l)
            self.l = l;
        end
        
        function self = set_d(d)
            self.d = d;
        end
        
        function self = set_u(u)
            self.u = u;
        end
    end
    
    methods (Static)
        % Initial concentration through cell, clean water input
        function y = cdf_ic(t, l, d, u)
            % Integral concentration left in the interval (-infinity, l]
            tmp_sqrt_d_t = sqrt(d .* t);
            tmp_l_sqrt_pi = (l .* sqrt(pi));
            tmp_tmp = 0.25 .* (l .^ 2 + u .^ 2 .* t .^ 2) ./ (d .* t);
            tmp_erf = 0.5 .* erf(0.5 .* (l - u .* t) ./ tmp_sqrt_d_t);
            y = 0.5 + ...
                tmp_sqrt_d_t .* exp(0.5 .* l .* u ./ d - tmp_tmp) ./ tmp_l_sqrt_pi - ...
                tmp_sqrt_d_t .* exp(0.25 .* l .^ 2 ./ (d .* t) - tmp_tmp) ./ tmp_l_sqrt_pi - ...
                u .* t .* tmp_erf ./ l - ...
                0.5 .* u .* t .* erf(0.5 .* u .* sqrt(t) ./ sqrt(d)) ./ l + ...
                tmp_erf;
            y(t == 0) = 1;
        end
        
        % Zero initial concentration through cell, instantaneous injection
        % of solute at t = 0
        function y = cdf_if(t, l, d, u)
            % Integral concentration left in the interval (-infinity, l]
            y = 0.5 .* (1 - ...
                2 .* sqrt(d) .* exp(0.5 .* l .* u ./ d - 0.25 .* (l .^ 2 + u .^ 2 .* t .^ 2) ./ (d .* t)) ./ (u .* sqrt(pi .* t)) + ...
                erf(0.5 .* (l - u .* t) ./ sqrt(d .* t)));
            y(t == 0) = 1;
        end
    
        function t_gen = rng(l, d, u, ic_type, varargin)
            if numel(varargin) == 0
                sz = size(l);
                num_cells = prod(sz);
                batch_size = 50;
                num_batches = ceil(num_cells / batch_size);
                ini_t = ones(sz);
                options = optimset('Display', 'off'); % , 'TolX', self.EPSILON
                t_gen = nan(sz);
                first_cell_idx = 1;
                rn = unifrnd(0, 1, sz);
                for idx = 1:num_batches
                    last_cell_idx = min(num_cells, first_cell_idx + batch_size - 1);
                    t_gen(first_cell_idx:last_cell_idx) = lsqnonlin(@f, ini_t(first_cell_idx:last_cell_idx), ...
                        zeros(last_cell_idx - first_cell_idx + 1, 1), inf(last_cell_idx - first_cell_idx + 1, 1) , options);
                    first_cell_idx = last_cell_idx + 1;
                end
            end
            
            return
            
            function y = f(t)
                l_int = l;
                d_int = d;
                u_int = u;
                if numel(l) > 1
                    l_int = l(first_cell_idx:last_cell_idx);
                end
                if numel(d) > 1
                    d_int = d(first_cell_idx:last_cell_idx);
                end
                if numel(u) > 1
                    u_int = u(first_cell_idx:last_cell_idx);
                end
                switch ic_type
                    case 1
                        y = solute_residence_generator.cdf_ic(t, l_int, d_int, u_int) - rn(first_cell_idx:last_cell_idx);
                    case 2;
                        y = solute_residence_generator.cdf_if(t, l_int, d_int, u_int) - rn(first_cell_idx:last_cell_idx);
                    otherwise
                        error('Error! Wrong IC type.');
                end
            end
        end
    end
end