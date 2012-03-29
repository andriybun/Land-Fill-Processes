function [spatial, is_landfill_array, column_height_array] = define_geometry()
    % Landfill size:
    x_len = 10;
    y_len = 15;
    z_depth = 10;

    % Structure defining landfill geometric properties:
    spatial = struct();
    spatial.dx = 1;
    spatial.dy = 1;
    spatial.dz = 1;
    spatial.xn = round(x_len / spatial.dx);
    spatial.yn = round(y_len / spatial.dy);
    spatial.zn = z_depth / spatial.dz;
    
    warning('Dealing with remainders, while computing numbers of columns, not yet implemented!');
    
    % 3D array defining cells belonging to a landfill body (1 - true, 0 -
    % false)
    is_landfill_array = ones(spatial.xn, spatial.yn, spatial.zn);
    
    %% STUB:
    is_landfill_array(1:4, 1:3, :) = 0;
    is_landfill_array(7:end, 8:end, 4:end) = 0;
    %% END STUB
    
    % 2D array defining column heights of size dx * dy in the landfill body
    % conducting liquid
    column_height_array = sum(is_landfill_array, 3) * spatial.dz;
    
    % Number of rectangular columns of size dx * dy in the landfill body
    % conducting liquid
    spatial.num_columns = nnz(column_height_array);
    
end