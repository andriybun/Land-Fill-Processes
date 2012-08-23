classdef water_particle_class

    %% Properties
    properties (Access = public)
        pos_arr;                % position of particle in the system, coordinates z, y, x (for 2-D case - z, x)
                                % array num_particles * ndims
        t_out_arr;              % time of particles' escapes
        active_arr;             % array denoting
        dv;                     % volume of a single particle
    end
    
    properties (Access = private)
        num_dims;               % dimensionality of the system
        num_columns;            % number of active columns
        num_particles;          % number of water particles considered at the moment
        init_coords;            % coords to be set for particles during initialization
        hydraulic_params;
        lognrnd_param_definer;
        
        % Log-normal parameters throughout the system
        mu;
        sigma;
        
        sz;
        zn;
        dz;
        t_prev;
    end
    
    %% Methods
    methods (Static)
        function result = sub_array(arr, idx)
            if size(idx, 2) == 2
                result = arr(sub2ind(size(arr), idx(:, 1), idx(:, 2)));
            elseif size(idx, 2) == 3
                result = arr(sub2ind(size(arr), idx(:, 1), idx(:, 2), idx(:, 3)));
            else
                result = [];
            end
        end
    end
    
    methods (Access = public)
        function self = water_particle_class(geom_params, dv, hydraulic_params, lognrnd_param_definer)
            self.dv = dv;
            self.num_dims = ndims(geom_params.is_landfill_array);
            self = self.set_init_coords(geom_params);
            self.hydraulic_params = hydraulic_params;
            self.lognrnd_param_definer = lognrnd_param_definer;
            
            self.sz = [geom_params.zn, geom_params.yn, geom_params.xn];
            self.zn = geom_params.zn;
            self.dz = geom_params.dz;
            self.t_prev = 0;
            self.num_particles = 0;
        end
        
        function [self, out_flux] = add_water_in(self, t, vol, se)
            % vol - volume of water per cell
            [self.mu, self.sigma] = self.lognrnd_param_definer.get_params(self.hydraulic_params.k_sat ./ self.dz, se);
            [self, out_flux] = self.update(t, se);
            
            % number of newly created particles per cell
            np = ceil(vol / self.dv);
            % total number of newly created 
            total_np = np * self.num_columns;
            self.num_particles = self.num_particles + total_np;
            
            % temporary arrays to be appended to properties of existing cells
            tmp_pos_arr = repmat(self.init_coords, [np, 1]);
            
            tmp_mu = self.sub_array(self.mu, tmp_pos_arr);
            tmp_sigma = self.sub_array(self.sigma, tmp_pos_arr);
            tmp_t_out_arr = t + lognrnd(tmp_mu, tmp_sigma, total_np, 1);
            
            % appending properties to the existing cells
            self.pos_arr = cat(1, self.pos_arr, tmp_pos_arr);
            self.t_out_arr = cat(1, self.t_out_arr, tmp_t_out_arr);
        end
        
        function [self, out_flux] = update(self, t, se)
            % Some indices of particles that change position
            cond_idx = (self.t_out_arr <= t);
            cond_idx_rep = repmat(cond_idx, [1, self.num_dims]);
            cond_idx_3d = logical(horzcat(cond_idx, zeros(size(cond_idx, 1), self.num_dims - 1)));
            
            % Number of particles that change position
            nel = sum(cond_idx);
            
            % Particles flow down, update position
            self.pos_arr(cond_idx_3d) = self.pos_arr(cond_idx_3d) + 1;
            if self.num_particles > 0
                part_out_idx = (self.pos_arr(:, 1) > self.zn);
                
                % Calculate out flux
                nel_out = sum(part_out_idx);
                out_flux = nel_out * self.dv;
                nel = nel - nel_out;
                self.num_particles = self.num_particles - nel_out;
                
                % Get rid of particles that left the system
                idx_out = find(part_out_idx);
                self.pos_arr(idx_out, :) = [];
                self.t_out_arr(part_out_idx, :) = [];
                cond_idx(part_out_idx) = [];
                cond_idx_rep(idx_out, :) = [];
                
                % Update exit time for particles that moved
                tmp_pos = reshape(self.pos_arr(cond_idx_rep), [], self.num_dims);
                
                tmp_mu = self.sub_array(self.mu, tmp_pos);
                tmp_sigma = self.sub_array(self.sigma, tmp_pos);
                self.t_out_arr(cond_idx) = t + lognrnd(tmp_mu, tmp_sigma, nel, 1);
            else
                out_flux = 0;
            end
            self.t_prev = t;
        end
        
        function conc = get_concentrations(self)
            conc = zeros(self.sz);
            for idx = 1:self.num_particles
                c_idx = self.pos_arr(idx, :);
                conc(c_idx) = conc(c_idx) + 1;
            end
            conc = conc * self.dv;
        end
    end
    
    methods (Access = private)
        function self = set_init_coords(self, geom_params)
            sz = size(geom_params.is_landfill_array);
            active_columns = sum(geom_params.is_landfill_array, 1);
            self.num_columns = nnz(active_columns);
            self.init_coords = zeros(self.num_columns, self.num_dims);
            
            if self.num_dims == 2
                sz = cat(2, sz, 1);
            end
            
            active_idx = 1;
            for j = 1:sz(2)
                for k = 1:sz(3)
                    if active_columns(1, j, k) > 0
                        if self.num_dims == 2
                            self.init_coords(active_idx, :) = [1, j];
                        else
                            self.init_coords(active_idx, :) = [1, j, k];
                        end
                        active_idx = active_idx + 1;
                    end
                end
            end
        end
        
    end
end


