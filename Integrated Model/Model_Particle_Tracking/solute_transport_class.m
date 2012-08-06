classdef solute_transport_class
    properties (Access = public)
        l;
        d;
        c_curr;              % initial concentration
    end
    
    properties (Access = private)
        domain_type;        % type of domain (1 - matrix, 2 - channel)
        ic_type;            % type of initial conditions (1 - initial concentration, or 2 - initial injection)
        is_initialized;     % 
        u_ini;
        u_prev;
        t_prev;
        t_delay;            % time delay if there was no flow at the beginning
        udt;                % auxiliary variable to store flow history
        pdf_handle;         % handle of PDF valid for our initial conditions
        domain_size;        % size of domain
        c_curr_prev;        % initial concentration on previous step
    end
    
    methods (Access = public)
        function self = solute_transport_class(l, d, domain_type, ic_type)
            self.l = l;
            self.d = d;
            self.c_curr = ones(size(l));             % initialize concentration
            self.c_curr_prev = self.c_curr;
            self.domain_size = size(l);
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
        end
        
        function c = get_c_curr(self)
            c = self.c_curr;
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
            
            t = t - self.t_delay(self.is_initialized);
            dt = t - self.t_prev;
            is_initial_step = (t == 0);
            
            self.c_curr = self.c_curr_prev;
            
            switch self.domain_type
                case 1
%                     result(self.is_initialized) = self.pdf_handle(t, self.u_ini(self.is_initialized));
                    self.c_curr = self.pdf_handle(t, self.u_ini(self.is_initialized));
                    result(is_initial_step) = 0;
                case 2
%                     if t == 0
%                         result(self.is_initialized) = self.pdf_handle(t, self.u_ini(self.is_initialized));
%                     else
%                         ratio = u ./ self.u_ini;
%                         % tweak time axis
%                         self.udt = self.udt + dt .* ratio;
%                         % adjust along y-axis along with changes in 
%                         result(self.is_initialized) = self.pdf_handle(self.udt, self.u_ini(self.is_initialized)) .* ratio;
%                     end
                    ratio = u ./ self.u_ini;
                    % tweak time axis
                    self.udt = self.udt + dt .* ratio;
                    self.c_curr = self.pdf_handle(self.udt, self.u_ini(self.is_initialized));
                    result(is_initial_step) = Inf;
            end
            
            
            result(~is_initial_step) = (self.c_curr(~is_initial_step) - self.c_curr_prev(~is_initial_step)) ./ dt(~is_initial_step);
%             result = max(-self.c_curr ./ dt, result);
%             
%             if t ~= 0
%                 if self.t_prev == 0
%                     self.c_curr = self.c_curr + 2 .* result .* dt;
%                 else
%                     self.c_curr = self.c_curr + result .* dt;
%                 end
%             end
            self.u_prev = u;
            self.t_prev = t;
            self.c_curr_prev = self.c_curr;
        end
    end
    
    methods (Access = private)        
        % Initial concentration through cell, clean water input
        function y = ini_concentration(self, t, u)
            l = self.l;
            d = self.d;
%             sqrt_t = sqrt(t);
%             sqrt_d = sqrt(d);
%             sqrt_pi = sqrt(pi);
%             tmp1 = (1 ./ 4) .* (l .^ 2 + u .^ 2 .* t .^ 2) ./ (t .* d);
%             y = -(1 ./ 2) .* (u .* erf((1 ./ 2) .* u .* sqrt_t ./ sqrt_d) + ...
%                 u .* erf((1 ./ 2) .* (l - u .* t) ./ (sqrt_t .* sqrt_d)) + ...
%                 sqrt_d .* exp((1 ./ 4) .* l .^ 2 ./ (t .* d) - tmp1) ./ (sqrt_t .* sqrt_pi) - ...
%                 sqrt_d .* exp((1 ./ 2) .* l .* u ./ d - tmp1) ./ (sqrt_t .* sqrt_pi));
%             y(t == 0) = Inf;
            % Integral concentration left in the interval (-infinity, l]
            y = 1 ./ 2 .* (1 + ...
                2 .* sqrt(d .* t) .* exp((1 ./ 2) .* l .* u ./ d - (1 ./ 4) .* (l .^ 2 + u .^ 2 .* t .^ 2) ./ (d .* t)) ./ (l .* sqrt(pi)) - ...
                2 .* sqrt(d .* t) .* exp((1 ./ 4) .* l .^ 2 ./ (d .* t) - (1 ./ 4) .* (l .^ 2 + u .^ 2 .* t .^ 2) ./ (d .* t)) ./ (l .* sqrt(pi)) - ...
                u .* t .* erf((1 ./ 2) .* (l - u .* t) ./ sqrt(d .* t)) ./ l - ...
                u .* t .* erf((1 ./ 2) .* u .* sqrt(t) ./ sqrt(d)) ./ l + ...
                erf((1 ./ 2) .* (l - u .* t) ./ sqrt(d .* t)));
            y(t == 0) = 1;
        end
        
        % Zero initial concentration through cell, instantaneous injection
        % of solute at t = 0
        function y = instant_influx(self, t, u)
            l = self.l;
            d = self.d;
%             y = (1 ./ 4) .* (-l .* u .* t + 2 .* d .* t - l .^ 2) .* exp(-(1 ./ 4) .* (l - u .* t) .^ 2 ./ (d .* t)) ./ ...
%                 (l .* u .* sqrt(pi .* d .* t .^ 5));
%             y(t == 0) = 0;
            % Integral concentration left in the interval (-infinity, l]
            y = 1 ./ 2 .* (1 - ...
                2 .* sqrt(d) .* exp((1 ./ 2) .* l .* u ./ d - (1 ./ 4) .* (l .^ 2 + u .^ 2 .* t .^ 2) ./ (d .* t)) ./ (u .* sqrt(pi .* t)) + ...
                erf((1 ./ 2) .* (l - u .* t) ./ sqrt(d .* t)));
            y(t == 0) = 1;
        end
    end
end