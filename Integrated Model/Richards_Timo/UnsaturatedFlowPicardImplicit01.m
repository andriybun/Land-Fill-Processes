%% Template for a Simulation Model
%
% Example template to illustrate an algorithmic approach to running simulations
% in Matlab
%
% Author: T.J. Heimovaara
% Date: April 16, 2012
% version: 0.01

function UnsaturatedFlowPicardImplicit01
    clear
    close all

    %% Discretization
    %Initialize Solver (discretization of space and time)

    % Spatial Discretization
    % x,y,z
    % Define internodal boundaries (include at least all soillayer boundaries)
    zbot = -20;
    zin = (0:-0.2:zbot)';  %soil profile until 5 meters depth
    nin = length(zin);
    zn(1,1) = zin(1,1);
    zn(2:nin-2,1) = (zin(2:nin-2)+zin(3:nin-1))./2;
    zn(nin-1,1) = zin(nin);
    nn = length(zn);

    dzn = zn(2:end,1)-zn(1:end-1,1);
    dzin = zin(2:end,1)-zin(1:end-1,1);

    %Store model dimensions in ModelDim
    ModelDim.zn=zn;
    ModelDim.zin=zin;
    ModelDim.dzn=dzn;
    ModelDim.dzin=dzin;
    ModelDim.nn=nn;
    ModelDim.nin=nin;

    % boundary parameters
    %BoundaryPar.Tavg = 273+10;  %Average yearly temperature (K)
    %BoundaryPar.Trange = 15; %Temperature Amplitude (Kelvin)
    %BoundaryPar.tmin = 46; %Offset time for minimum temperature (from new year in day)
    BoundaryPar.qTop = @qBoundary;
    BoundaryPar.Ksurf = 1; % 1e-2;
    BoundaryPar.hamb = -0.95;

    %% Model initialization

    % model parameters
    %Nodal parameters
    alln = 1:nn;
    allin = 1:nin;

    hcap = 0.5; %Air entry head
    SoilPar.alpha(alln,1) = 1./hcap;
    SoilPar.n(alln,1) = 1.9;
    SoilPar.thetaR(alln,1) = 0.04;
    SoilPar.thetaS(alln,1) = 0.4;

    %This is a simple choice, however there is a relationship between nodes and
    %internodes!!!!
    SoilPar.alpha_in(allin,1) = 1./hcap;
    SoilPar.n_in(allin,1) = 1.9;
    SoilPar.thetaR_in(allin,1) = 0.04;
    SoilPar.thetaS_in(allin,1) = 0.4;
    SoilPar.Ksat_in(allin,1) = 1; %1e-2; %m/d

    % model initial states;
    rhow = 1000; % kg/m3 density of water
    g = 9.81; % m/s^2 gravitational acceleration
    zref = zbot; % Height of phreatic surface for initial condition
    hbot = BoundaryPar.hamb;
    hIni(alln,1)= hbot-(zn-zref);

    %% Dynamic loop (non-stationary model) %% Implicit solution with ode15i
    % Time Discretization
    % start time
    dtIni = 5;
    trange = [0:dtIni:3000]; % [0:1:4 5:5:550];
    dtmax = 20;

    tic

    %Initialize the statevariables for the simulation
    t = trange(1);
    hSim = hIni;
    States = Waterstates(t,hSim,SoilPar,ModelDim,BoundaryPar);
    theta = States.theta;
    q = States.q;
    Cm = Cm_function(hSim,SoilPar);

    niter = 0;
    maxiter = 45;
    miniter = 15;
    dt = dtmax;
    dtiter = dtmax;

    convcrit = 1e-7;

    %Initialize output matrices
    Time = t;
    hOut = hSim';
    qOut = q';
    thetaOut = theta';
    CumMBal = 0;
    CumMBalOut = CumMBal;
    dttOut = dt;

    while abs(t-trange(end)) > eps,
        %our primary variable will be updated during the iteration. Store the
        %current value for the previous time
        hSim_o = hSim;
        theta_o = theta;

        %initialize local iteration counter (for reducing timestep)
        niter = 0;
        converged = false;

        %Check timestep
        %output time
        dtoutr=(t-trange);
        dtoutidx = find(dtoutr<0);
        next_tout = trange(dtoutidx(1));
        dtout = next_tout-t;

        %max time step related to flow rate
        % 50% change in pressure head
        change = 0.9;
        dtflow1 = min(abs((change.*Cm*hSim.*dzin)./-(q(2:nin,1)-q(1:nin-1,1))));
        dtflow2 = min(abs((change.*theta_o.*dzin)./-(q(2:nin,1)-q(1:nin-1,1))));

        dt = min([dtflow2,dtmax,dtiter,dtout]);
        %dt = min([dtmax,dtiter,dtout]);

        while ~converged,
            niter = niter + 1;
            Cm = Cm_function(hSim,SoilPar);
            Km = Kmat(hSim,SoilPar,BoundaryPar,ModelDim);
            Yv = Yvec(t,hSim,SoilPar,BoundaryPar,ModelDim);

            %(Km+Cm)*delta = Km*h+Yv-(theta-theta_o)./dt
            %Ax =b --> x = A\b

            Atmp = (-Km+Cm./dt);
            btmp = Km*hSim + Yv - (theta-theta_o)./dt;

            delta = Atmp\btmp;

            %disp([niter t dt log10(max(abs(delta)))])

            hSim = hSim + delta;
            States = Waterstates(t,hSim,SoilPar,ModelDim,BoundaryPar);
            theta = States.theta;
            q = States.q;

            if max(abs(delta))<convcrit,
                converged = true;
                if niter >= maxiter,
                    dtiter = dtiter./3;
                elseif niter <= miniter,
                    dtiter = dtiter.*2;
                end
                t = t+dt;
                %Calculate Cumulative MassBalance
                TotWatMass_o = -sum(theta_o(:,1).*dzin);
                TotWatMass = -sum(theta(:,1).*dzin); 
                NFlux = -(q(1,1)-q(end,1));
                % Total change of mass in column should be equal to total the difference
                % between the time integrated flux at the top and the bottom boundaries

                MassBal = (TotWatMass-TotWatMass_o)-NFlux.*dt;
                CumMBal = CumMBal+MassBal;

                %%Check for output
                if abs(t-next_tout)<eps,
                    Time = [Time;t];
                    hOut = [hOut;hSim'];
                    qOut = [qOut;q'];
                    thetaOut = [thetaOut;theta'];
                    CumMBalOut = [CumMBalOut; CumMBal];
                    dttOut = [dttOut;dt];
                end

            end
        end
    end

    toc

    minX = 41;
    step = 20;
    maxX = 101;
    factor = sum(qOut(:, 1)) * dtIni;
    if factor == 0
        factor = 1;
    end
    out_flux_normalized = qOut(:, minX:step:maxX) / factor;
    plot(Time, out_flux_normalized);
    sum(qOut(:, minX:step:maxX))'

    [mu, sigma] = pick_params(Time, out_flux_normalized);
    dif_prev = inf;
    for idx = 1:numel(mu)
        iters = 0;
        while true && (iters < 1000)
            sigma(idx) = sigma(idx) + 2e-3;
            out_flux_lognorm = out_flux_lognrnd_pdf(Time, mu(idx), sigma(idx));
            out_flux_lognorm(1) = 0;
            dif = square_diff(out_flux_normalized(:, idx), out_flux_lognorm);
            if dif > dif_prev
                break
            end
            dif_prev = dif;
            iters = iters + 1;
        end
    end
    hold on;
    plot(Time, out_flux_lognrnd_pdf(Time, mu, sigma), 'Color', [1, 0.5, 0.2]);
    hold off;
    
    figure(2);
    subplot(2,1,1);
    plot(mu);
    subplot(2,1,2);
    plot(sigma);
    
    figure(1)
    clf
    hold on
    plot(hOut',zn);
    %plot temperature as a function of depth for selected times
    figure(2)
    clf
    hold on
    plot(Time,hOut)
    
    return

    function [mux, sigmax] = pick_params(t, out_flx)
        t = repmat(t, [1, size(out_flx, 2)]);
        dt = t(2:end, :) - t(1:end-1, :);
        mn = sum(t(2:end, :) .* (out_flx(1:end-1, :) + out_flx(2:end, :)) ./ 2 .* dt, 1);
        mnr = repmat(mn, [size(out_flx, 1) - 1, 1]);
        varn = sum((t(2:end, :) - mnr).^2 .* (out_flx(1:end-1, :) + out_flx(2:end, :)) ./ 2 .* dt);
        params = get_mu_sigma(cat(1, mn, varn));
        mux = params(1, :);
        sigmax = params(2, :);
    end

    function params = get_mu_sigma(moments)
        m = squeeze(moments(1, :));
        v = squeeze(moments(2, :));
        sigmax = sqrt(log(v ./ (m .* m) + 1));
        mux = log(m) - sigmax .* sigmax / 2;
        params = cat(1, mux, sigmax);
    end

    function res = out_flux_lognrnd_pdf(t, muv, sigmav)
        num = numel(muv);
        res = zeros(numel(t), num);
        for idx = 1:num
            res(:, idx) = 1 ./ (sqrt(2 * pi) .* sigmav(idx) .* t) .* exp(-(log(t) - muv(idx)).^2 ./ (2 .* sigmav(idx).^2));
            res(1, idx) = 0;
        end
    end

    function [res, ratio] = square_diff(v1, v2, ratiox)
        if nargin < 3
            mv1 = max(v1);
            mv2 = max(v2);
            ratio = mv1 / mv2;
        else
            ratio = ratiox;
        end
        res = sum((v1 - v2 .* ratio) .* (v1- v2 .* ratio));
    end
end