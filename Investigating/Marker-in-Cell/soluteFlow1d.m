function qH = soluteFlow1d(t, C, SoilPar, ModelDim, BoundaryPar)
    % Function for calculating rates for a Heat Flow in a spatial domain
    % qH are the heat fluxes at all internodal points in the spatial domain;
    % t is time
    % T are the temperatures in all cells in the spatial domain
    % SoilPar contains the Model Parameters (thermal diffusivity of all
    % InterNodal points)
    % ModelDim contains the ModelDiscretization (coordinates of nodes and 
    % internodes, deltas of nodal points and
    % internodal points)

    %Copy parameter to local variable
    nn = ModelDim.nn;
    nin = ModelDim.nin;
    dzn = ModelDim.dzn;
    dzin = ModelDim.dzin;

    %Diffusive Flux
    D = SoilPar.D(1:ModelDim.nin, 1);               % [m²/s]
    thetain = SoilPar.thetain(1:ModelDim.nin, 1);     % [m³/m³]
    
    %Concentration boundary conditions upper/left boundary condition
    C(1) = BoundaryPar.CTop;

    qD(1, 1) = -D(1, 1) .* thetain(1, 1) .* (C(1, 1) ./ dzn(1));
    qD(2:nin-1, 1) = -D(2:nin-1, 1) .* thetain(2:nin-1, 1) .* (C(2:nn, 1) - C(1:nn-1, 1)) ./ dzn(1:nn-1);

    % Lower/Right boundary conditions
    % zero concentration gradient ...
    qD(nin,1)= 0;

    qH = qD;
end