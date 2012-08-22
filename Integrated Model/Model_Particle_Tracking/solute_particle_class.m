classdef solute_particle_class

    %% Properties
    properties (Access = public)
        state_arr;              % array of state of particle (1 - absorbed(default), 2 - released, 3 - free)
                                % depending on state, different initial conditions are applied:
                                % (0 - none (the particle stays), 1 - constant initial concentration, 
                                % 2 - instantaneous influx)
        pos_arr;                % position of particle in the system, coordinates z, y, x (for 2-D case - z, x)
                                % array num_particles * ndims
        t_out_arr;              % time of particles' escapes
        l_arr;                  % height of cells in which particles reside
        active_arr;             % array denoting
        d;
    end
    
    properties (Access = private)
        num_particles;
        ddims;
        dsize;
        residence_time_generator;
        l;
        coords;
        t_prev;
        % State transition params - probabilities that particles change
        % state in unit time
        prob_abs_rel = 5e-2;
        prob_free_abs = 3e-2;
        lambda_abs_rel;
        lambda_free_abs;
        state_transition;
    end
    
    %% Methods
    
    %% Public methods
    methods (Access = public)
        
        % Constructor
        function self = solute_particle_class(geom_params, solute_params)
            num_particles_arr = ceil(geom_params.c_ini / solute_params.dv);
            self.num_particles = sum(num_particles_arr);
            sz = size(num_particles_arr);
            self.ddims = ndims(num_particles_arr);
            self.dsize = sz;
            
            self.l = geom_params.l;
            self.coords = cumsum(geom_params.l);
            self.d = solute_params.d;
            self.t_prev = 0;
            
            % Initialize properties of solute particles
            self.state_arr = ones(self.num_particles, 1);
            self.pos_arr = zeros(self.num_particles, ndims(num_particles_arr));
            self.t_out_arr = inf(self.num_particles, 1);
            self.l_arr = zeros(self.num_particles, 1);
            self.active_arr = true(self.num_particles, 1);
            
            part_idx = 1;                       % index of first particle of batch for current cell
            switch ndims(num_particles_arr)
                % 2-D case
                case 2
                    for i = 1:sz(1)
                        for j = 1:sz(2)
                            batch_size = num_particles_arr(i, j);
                            self.pos_arr(part_idx:part_idx + batch_size - 1, :) = repmat([i, j], [batch_size, 1]);
                            self.l_arr(part_idx:part_idx + batch_size - 1) = self.l(i, j);
                            part_idx = part_idx + batch_size;
                        end
                    end
                    % 3-D case
                case 3
                    for i = 1:sz(1)
                        for j = 1:sz(2)
                            for k = 1:sz(3)
                                batch_size = num_particles_arr(i, j, k);
                                self.pos_arr(part_idx:part_idx + batch_size) = repmat([i, j, k], [batch_size, 1]);
                                self.l_arr(part_idx:part_idx + batch_size - 1) = self.l(i, j, k);
                                part_idx = part_idx + batch_size;
                            end
                        end
                    end
                otherwise
                    error('Error: Invalid dimensions.');
            end
            
            % State switching params
            self.lambda_abs_rel = -log(1 - self.prob_abs_rel);
            self.lambda_free_abs = -log(1 - self.prob_free_abs);
        end
        
        function [self, num_escaped] = update(self, t, u)
            if numel(u) == 1
                u = repmat(u, size(self.l_arr));
            elseif numel(u) == numel(self.l)
                u_in = u;
                u = zeros(size(self.l_arr));
                tmp_idx_all = self.pos_arr(:, 1);
                u(self.active_arr) = u_in(tmp_idx_all(self.active_arr));
            end
            if numel(self.d) == 1
                d_in = repmat(self.d, size(self.l_arr));
            else
                d_in = self.d;
            end
            self = self.set_state_transition_matrix(t);
            % If current time is after t_out_arr we change location of particles, states
            cond_idx = self.active_arr & (self.t_out_arr <= t);
            cond_idx_3d = logical(horzcat(cond_idx, zeros(size(cond_idx, 1), self.ddims - 1)));
            % Particles flow down
            self.pos_arr(cond_idx_3d) = self.pos_arr(cond_idx_3d) + 1;
            self.active_arr(self.pos_arr > self.coords(end)) = false;
            % Indices of particles that left system on the last step
            just_escaped_idx = ~self.active_arr & logical(cond_idx);
            num_escaped = sum(just_escaped_idx);
            % Their state becomes free
            self.state_arr(cond_idx) = 3;
            % Generate new states of other particles
            self = self.switch_state(~cond_idx);
            % Generate new exit times for particles
            if self.ddims == 2
                self.l_arr(self.active_arr) = self.l(sub2ind(self.dsize, self.pos_arr(self.active_arr, 1), ...
                    self.pos_arr(self.active_arr, 2)));
            else
                self.l_arr(self.active_arr) = self.l(sub2ind(self.dsize, self.pos_arr(self.active_arr, 1), ...
                    self.pos_arr(self.active_arr, 2), self.pos_arr(self.active_arr, 3)));
            end
            self.t_out_arr((self.state_arr == 1) & (self.active_arr)) = inf;
            is_released_idx = (self.state_arr == 2) & (self.active_arr);
            self.t_out_arr(is_released_idx) = t + solute_residence_generator.rng(self.l_arr(is_released_idx), ...
                d_in(is_released_idx), u(is_released_idx), 1);
            is_free_idx = (self.state_arr == 3) & (self.active_arr);
            self.t_out_arr(is_free_idx) = t + solute_residence_generator.rng(self.l_arr(is_free_idx), ...
                d_in(is_free_idx), u(is_free_idx), 2);
            self.t_prev = t;
        end
        
    end
    
    %% Private methods
    methods (Access = private)
        function ic_type = get_ic_type(self)
            ic_type = nan(size(self.state_arr));
            ic_type(self.state_arr == 1) = 0;
            ic_type(self.state_arr == 2) = 1;
            ic_type(self.state_arr == 3) = 2;
        end
        
        function self = switch_state(self, idx)
            if nargin == 1
                self.state_arr = self.non_unif_rnd(self.state_transition(self.state_arr, :));
            else
                self.state_arr(idx) = self.non_unif_rnd(self.state_transition(self.state_arr(idx), :));
            end
        end
        
        function self = set_state_transition_matrix(self, t)
            dt = t - self.t_prev;
            % State transition
            a = 1 - exp(-self.lambda_abs_rel * dt);
            b = 1 - exp(-self.lambda_free_abs * dt);
            self.state_transition = [ ...
                1 - a       a           0;
                0           (dt == 0) 	(dt ~= 0);
                b           0           1 - b ];
        end
    end
       
    methods (Static)
        % Abstract stuff
        function rn = non_unif_rnd(prob_matr)
            % Preprocess input parameters
            if find(size(prob_matr) == 1) > 0
                % If probability array is vector, make it row vector
                prob_matr = reshape(prob_matr, [1, numel(prob_matr)]);
            end
            % Generate random number for each row
            sz = [size(prob_matr, 1), 1];
            
            % Initialize resulting matrix
            rn = nan(sz);
            
            % Generate random numbers
            cdf = horzcat(zeros(sz), cumsum(prob_matr, 2));
            if find(cdf(:, end) ~= 1) > 1
                error('Error: Wrong probability distribution. Sum of all probabilities is %3.3f instead of 1\n', cdf(end));
            end
            
            tmp_rn = rand(sz);
            for idx = 1:size(cdf, 2)-1
                cond_idx = (tmp_rn > cdf(:, idx));
                rn(cond_idx) = idx;
            end
        end
    end
end