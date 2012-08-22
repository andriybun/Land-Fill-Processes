classdef solute_transport_class
    properties (Access = public)
        l;
        d;
        c_ini;              % initial concentration
        c_curr;             % current concentration (relative to the initial)
    end
    
    properties (Access = private)
        EPSILON;
        domain_type;        % type of domain (1 - matrix, 2 - channel)
        ic_type;            % type of initial conditions (1 - initial concentration, or 2 - initial injection)
        is_active;
        is_initialized;     % 
        u_ini;
        u_prev;
        t_prev;
        t_delay;            % time delay if there was no flow at the beginning
        udt;                % auxiliary variable to store flow history
        pdf_handle;         % handle of PDF valid for our initial conditions
        domain_size;        % size of domain
        c_prev;             % concentration on previous step
    end
    
    methods (Access = public)
        function self = solute_transport_class(l, d, c_ini, domain_type, ic_type)
            self.EPSILON = 1e-7;
            if nargin > 0
                self.l = l;
                if numel(d) == 1
                    self.d = repmat(d, size(l));
                else
                    self.d = d;
                end
                self.c_ini = c_ini;
                self.c_curr = ones(size(l));             % initialize concentration
                self.c_prev = self.c_curr;
                self.domain_size = size(l);
                self.is_active = true(self.domain_size);
                self.is_initialized = false(self.domain_size);
                self.udt = zeros(self.domain_size);
                self.u_ini = nan(self.domain_size);
                self.ic_type = ic_type;
                self.domain_type = domain_type;
                self.t_delay = nan(self.domain_size);
                self.t_prev = 0;
                switch ic_type
                    case 1
                        self.pdf_handle = @self.ini_concentration;
                    case 2;
                        self.pdf_handle = @self.instant_influx;
                    otherwise
                        error('Error! Wrong IC type.');
                end
            else
                self.is_initialized = false;
            end
        end
        
        function self = set_c_ini(self, c_ini)
            self.c_ini = c_ini;
        end
        
        function c = get_c_ini(self)
            c = self.c_ini;
        end
        
        function c = get_c_curr(self)
            c = self.c_curr .* self.c_ini;
        end
        
        function [self, result] = flush(self, t, u)
            result = zeros(size(u));
            
            if numel(self.is_initialized) == 1 && (~self.is_initialized)
                self.u_ini = zeros(self.domain_size);
                self.u_prev = zeros(self.domain_size);
                self.is_initialized = (u > 0);
                self.t_delay = t .* ones(self.domain_size);
                self.u_ini(self.is_initialized) = u(self.is_initialized);
            else
                idx_initialize = logical(~self.is_initialized .* (u > 0));
                self.u_ini(idx_initialize) = u(idx_initialize);
                self.is_initialized(idx_initialize) = true;
                self.t_delay(idx_initialize) = t;
            end
            
            t = repmat(t, self.domain_size);
            t(self.is_initialized) = t(self.is_initialized) - self.t_delay(self.is_initialized);
            dt = t - self.t_prev;
            is_initial_step = (t == 0);
            
            self.c_curr = self.c_prev;
            
            comp_idx = self.is_initialized & self.is_active;
            
            switch self.domain_type
                case 1
                    self.c_curr(comp_idx) = self.pdf_handle(t(comp_idx), self.u_ini(comp_idx), comp_idx);
                    result(is_initial_step) = 0;
                case 2
                    ratio = nan(self.domain_size);
                    ratio(comp_idx) = u(comp_idx) ./ ...
                        self.u_ini(comp_idx);
                    % tweak time axis
                    self.udt(comp_idx) = self.udt(comp_idx) + dt(comp_idx) .* ratio(comp_idx);
                    self.c_curr(comp_idx) = self.pdf_handle(self.udt(comp_idx), self.u_ini(comp_idx), comp_idx);
            end
            
            result(~is_initial_step) = (self.c_curr(~is_initial_step) - self.c_prev(~is_initial_step));
            self.is_active(self.c_curr < self.EPSILON) = false;
            
            result = result .* self.c_ini;
            
            self.u_prev = u;
            self.t_prev = t;
            self.c_prev = self.c_curr;
        end
    end
    
    methods (Access = private)        
        % Initial concentration through cell, clean water input
        function y = ini_concentration(self, t, u, idx)
            if nargin < 4
                l = self.l;
                d = self.d;
            else
                l = self.l(idx);
                d = self.d(idx);
            end
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
        function y = instant_influx(self, t, u, idx)
            if nargin < 4
                l = self.l;
                d = self.d;
            else
                l = self.l(idx);
                d = self.d(idx);
            end
            % Integral concentration left in the interval (-infinity, l]
            y = 0.5 .* (1 - ...
                2 .* sqrt(d) .* exp(0.5 .* l .* u ./ d - 0.25 .* (l .^ 2 + u .^ 2 .* t .^ 2) ./ (d .* t)) ./ (u .* sqrt(pi .* t)) + ...
                erf(0.5 .* (l - u .* t) ./ sqrt(d .* t)));
            y(t == 0) = 1;
        end
    end
end