function velocity_array = generate_velocities(spatial, column_height_array, mu, sigma)
    velocity_array = zeros(size(column_height_array));
    velocity_array(column_height_array > 0) = lognrnd(mu, sigma, spatial.num_columns, 1);
end