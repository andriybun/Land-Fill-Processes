function nf = NettoFlux(t, c, SoilPar, ModelDim, BoundaryPar)
    qH = SoluteFlow1d(t, c, SoilPar, ModelDim, BoundaryPar);
    nin = ModelDim.znin;
    dzin = ModelDim.dzin;

    nf = -(qH(2:nin) - qH(1:nin-1)) ./ (dzin(1:nin-1));
    
    return
    
    function qD = SoluteFlow1d(t, c, SoilPar, ModelDim, BoundaryPar)
        nn = ModelDim.znn;
        nin = ModelDim.znin;
        dzn = ModelDim.dzn;
        
        % Diffusive Flux
        d = SoilPar.d(1:ModelDim.znin, 1);                  % [m²/s]
        thetaIn = SoilPar.thetaIn(1:ModelDim.znin, 1);      % [m³/m³]
        
        % Concentration boundary conditions upper/left boundary condition
        c(1) = BoundaryPar.cTop;
        
        % Internodal flux
        qD = nan(nin, 1);
        qD(1) = -d(1) .* thetaIn(1) .* (c(1) ./ dzn(1));
        qD(2:nin-1) = -d(2:nin-1) .* thetaIn(2:nin-1) .* (c(2:nn) - c(1:nn-1)) ./ dzn(1:nn-1);
        qD(nin) = 0;
    end
end