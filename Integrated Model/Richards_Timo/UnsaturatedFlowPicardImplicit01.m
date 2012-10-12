%% Template for a Simulation Model
%
% Example template to illustrate an algorithmic approach to running simulations
% in Matlab
%
% Author: T.J. Heimovaara
% Date: April 16, 2012
% version: 0.01
clear
close all

%% Discretization
%Initialize Solver (discretization of space and time)

% Spatial Discretization
% x,y,z
% Define internodal boundaries (include at least all soillayer boundaries)
zin = (0:-0.01:-1)';  %soil profile until 5 meters depth
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
BoundaryPar.qTop = @qBoundary; %.
BoundaryPar.Ksurf = 1e-2;
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
SoilPar.Ksat_in(allin,1) = 1e-2; %m/d

% model initial states;
rhow = 1000; % kg/m3 density of water
g = 9.81; % m/s^2 gravitational acceleration
zref = -1.0; % Height of phreatic surface for initial condition
hbot = -2.5;
hIni(alln,1)= hbot-(zn-zref);

%% Dynamic loop (non-stationary model) %% Implicit solution with ode15i
% Time Discretization
% start time
trange = [0:1:4 5:5:250];%0:5:250*365; %5 years)
dtmax = 10;

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

figure(1)
clf
hold on
plot(hOut',zn);
%plot temperature as a function of depth for selected times
figure(2)
clf
hold on
plot(Time,hOut)

