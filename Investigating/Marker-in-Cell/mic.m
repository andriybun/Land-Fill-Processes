function [t, cMark, cNode, zMark, nMarkOut] = mic(t, nOut, zMark, cMark,  v, trange, modelPar, geometry, boundaryPar)
    % Marker-in-cell function

    %% TMP workaround
    if ~isfield(geometry, 'z')
        geometry.z = geometry.zn;
    end
    %% end

    % Get velocity at markers' coordinates
    vMark = interp1(geometry.z, v, zMark, 'pchip', 'extrap');
    
    % Max timestep is based on cell size and flow velocity
    dtMax = 0.2 .* max(abs(geometry.dzn)) ./ max(abs(v));
    
    % Advect markers
    zMark = zMark + vMark .* dtMax;
    
    % Markers leaving system
    idxOut = (zMark < geometry.z(end));
    nMarkOut = sum(idxOut);
    zMark(idxOut) = zMark(idxOut) - geometry.z(end);
    cMark(idxOut) = boundaryPar.CTop;

    % Estimate new nodal concentrations by interpolation
    [cAvg, sIdx] = calcNodalValues(zMark, cMark, geometry);

    % Rates of change
    dCdt = nettoFlux(t, cAvg, modelPar, geometry, boundaryPar);

    % Time step should be small enough to prevent negative states from
    % occuring and should allow for all output times to be reached
    % Max time step depends on flow rate and on rates of change due to
    % processes
    dttest = max(abs(0.5 .* cAvg ./ dCdt)); %criterium still needs to be implemented
    dtout = trange(nOut + 1) - t;
    delt = min([dttest(:)' dtMax dtout]);

    % RK4 steps for non-advective transport processes
    k1 = delt .* dCdt;
    k2 = delt .* nettoFlux(t+delt/2, cAvg+k1/2, modelPar, geometry, boundaryPar);
    k3 = delt .* nettoFlux(t+delt/2, cAvg+k2/2, modelPar, geometry, boundaryPar);
    k4 = delt .* nettoFlux(t+delt, cAvg+k3, modelPar, geometry, boundaryPar);

    % Update states
    % Advection and Diffusion
    cNode = cAvg + (k1 + 2 * k2 + 2 * k3 + k4) / 6;

    % Calculate sub-grid diffusion term (Gerya 2010, chapter 10).
    % interpolate new temperature to Markers
    cmTmp(sIdx, 1) = interp1(geometry.zn, cNode, zMark(sIdx), 'pchip', 'extrap');

    d = 0.8;
    dm = interp1(geometry.zin, modelPar.D, zMark(sIdx), 'pchip', 'extrap');

    dzm = max(abs(geometry.dzn));
    dtDiff = 1 ./ (dm .* 2 .* dzm.^2);

    dCsgm = (cmTmp - cMark) .* (1 - exp(-d .* delt ./ dtDiff));
    dCsubgrid = interp1(zMark(sIdx), dCsgm, geometry.zn, 'pchip', 'extrap');

    dCsg2 = interp1(geometry.zn, dCsubgrid, zMark(sIdx), 'pchip', 'extrap');
    cMark = cMark + dCsg2;

    %update time
    t = t + delt;

end