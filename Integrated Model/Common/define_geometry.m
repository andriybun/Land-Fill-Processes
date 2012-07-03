function spatial_params = define_geometry()
    % Landfill size:
    x_len = 1;%10;
    y_len = 1;%15;
    z_depth = 10;

    % Structure defining landfill geometric properties:
    spatial_params = struct();
    spatial_params.dx = 1;
    spatial_params.dy = 1;
    spatial_params.dz = 1;
    spatial_params.xn = round(x_len / spatial_params.dx);
    spatial_params.yn = round(y_len / spatial_params.dy);
    spatial_params.zn = z_depth / spatial_params.dz;
    
    warning('Dealing with remainders, while computing numbers of columns, not yet implemented!');
    
    % 3D array defining cells belonging to a landfill body (1 - true, 0 -
    % false)
    spatial_params.is_landfill_array = ones(spatial_params.xn, spatial_params.yn, spatial_params.zn);
    
%     %% STUB:
%     spatial_params.is_landfill_array(1:4, 1:3, :) = 0;
%     spatial_params.is_landfill_array(7:end, 8:end, 4:end) = 0;
%     %% END STUB
    
    % 2D array defining column heights of size dx * dy in the landfill body
    % conducting liquid
    spatial_params.column_height_array = sum(spatial_params.is_landfill_array, 3) * spatial_params.dz;
    
    % Number of rectangular columns of size dx * dy in the landfill body
    % conducting liquid
    spatial_params.num_columns = nnz(spatial_params.column_height_array);
    % Area of top of the landfill
    spatial_params.area = spatial_params.num_columns * spatial_params.dx * spatial_params.dy;
end