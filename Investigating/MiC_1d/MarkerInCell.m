function [t, cMark, cNode, zMark, nMarkOut, massOut] = MarkerInCell(t, iTime, zMark, cMark,  v, tRange, SoilPar, ModelDim, BoundaryPar)
    %
    % Marker-in-cell function
    %
    
    % EPSILON = 1e-7;
    
    % Max timestep is based on cell size and flow velocity
    dtMax = max(abs(ModelDim.dzn)) ./ max(abs(v));
    dtOut = tRange(iTime + 1) - t;
    deltaT = min(dtMax, dtOut);
    
    % Estimate new nodal concentrations by interpolation
%     [cAvg, ~] = ComputeNodalValues(zMark, cMark, ModelDim);
    [cAvg, ~] = ComputeNodalValues([ModelDim.zin(1); zMark], [BoundaryPar.cTop; cMark], ModelDim);

    % Change of concentration due to non-advective transport processes
    % Numerically solve diffusion equation. Using Runge-Kutta 4th order
    % scheme:
    k1 = deltaT .* NettoFlux(t, cAvg, SoilPar, ModelDim, BoundaryPar);
    k2 = deltaT .* NettoFlux(t + deltaT/2, cAvg + k1/2, SoilPar, ModelDim, BoundaryPar);
    k3 = deltaT .* NettoFlux(t + deltaT/2, cAvg + k2/2, SoilPar, ModelDim, BoundaryPar);
    k4 = deltaT .* NettoFlux(t + deltaT, cAvg + k3, SoilPar, ModelDim, BoundaryPar);
    cNode = cAvg + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    
%     cNode(1) = BoundaryPar.cTop;
    
    % Calculate sub-grid diffusion term (Gerya 2010, chapter 10).
    % interpolate new temperature to Markers
    dC = cNode - cAvg;
    cAvgMark = ComputeMarkerValues(ModelDim.zn, cAvg, zMark);

    d = 1;
    dMark = ComputeMarkerValues(ModelDim.zin, SoilPar.d, zMark);

    dzMax = max(abs(ModelDim.dzin));
    dtDiff = 1 ./ (dMark .* 2 ./ dzMax.^2);

    dcSubgridMark = (cAvgMark - cMark) .* (1 - exp(-d .* deltaT ./ dtDiff));
    dcSubgridNode = ComputeNodalValues(zMark, dcSubgridMark, ModelDim);
    dcRemaining = dC - dcSubgridNode;
    dcRemainingMark = ComputeMarkerValues(ModelDim.zn, dcRemaining, zMark);
    cMark = cMark + dcSubgridMark + dcRemainingMark;
    
    % Get velocity at markers' coordinates
    vMark = ComputeMarkerValues(ModelDim.zn, v, zMark);
    
    dzMark = vMark .* deltaT;
    zMarkNew = zMark + dzMark;
    idxOut = (zMarkNew < ModelDim.zn(end));
    nMarkOut = sum(idxOut);
    nMarkPerNode = numel(cMark) / ModelDim.znn;
    massOut = 0;
    
    if any(idxOut)
        dtOut = (ModelDim.zn(end) - zMark(idxOut)) ./ vMark(idxOut);
        dtRemain = deltaT - dtOut;
        zMarkNew(idxOut) = ModelDim.zin(1) + v(1) * dtRemain;
        massOut = mean(cMark(idxOut)) * abs(ModelDim.dzin(ModelDim.znin - 1)) * nMarkOut / nMarkPerNode;
        cMark(idxOut) = BoundaryPar.cTop;
    end
    zMark = zMarkNew;
    
    cNode = ComputeNodalValues(zMark, cMark, ModelDim);
    
    % Update time
    t = t + deltaT;
    
    return
    
    function MdOut = Expand(Md)
        MdOut.znn = Md.znn + 1;
        MdOut.znin = Md.znin + 1;
        MdOut.dzin = cat(1, Md.dzin(1), Md.dzin);
        MdOut.zin = cat(1, Md.zin(1) - MdOut.dzin(1), Md.zin);
        MdOut.dzn = cat(1, Md.dzin(1), Md.dzn);
        MdOut.zn = cat(1, Md.zn(1) - MdOut.dzn(1), Md.zn);
    end
end