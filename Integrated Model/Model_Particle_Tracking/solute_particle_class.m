classdef solute_particle_class

    %% Properties
    properties (Access = public)
        state_arr;              % array of state of particle (1 - absorbed(default), 2 - released, 3 - free)
                                % depending on state, different initial conditions are applied:
                                % (0 - none (the particle stays), 1 - constant initial concentration, 
                                % 2 - instantaneous influx)
        pos_arr;                % position of particle in the system, coordinates x, y, z (for 2-D case - x, z)
                                % array num_particles * ndims
        t_out_arr;              % time of particles' escapes
        l_arr;                  % height of cells in which particles reside
    end
    
    properties (Access = private)
        num_particles;
        siz;
        residence_time_generator;
        l;
        % State transition params
        prob_abs_rel = 1e-1;
        prob_free_abs = 5e-2;;
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
            self.siz = sz;
            
            self.l = geom_params.l;
            
            % Initialize properties of solute particles
            self.state_arr = ones(self.num_particles, 1);
            self.pos_arr = zeros(self.num_particles, ndims(num_particles_arr));
            self.t_out_arr = inf(self.num_particles, 1);
            self.l_arr = zeros(self.num_particles, 1);
            
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
            
            % State transition
            self.state_transition = [ ...
                1 - self.prob_abs_rel   self.prob_abs_rel   0;
                0                       0                   1;
                self.prob_free_abs      0                   1 - self.prob_free_abs ];
        end
        
        function self = switch_state(self)
            self.state_arr = self.non_unif_rnd(self.state_transition(self.state_arr, :));
        end
    end
    
    %% Private methods
    methods (Access = private)
        function ic_type = get_ic_type(self, state)
            ic_type = nan(size(state));
            ic_type(state == 1) = 0;
            ic_type(state == 2) = 1;
            ic_type(state == 3) = 2;
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