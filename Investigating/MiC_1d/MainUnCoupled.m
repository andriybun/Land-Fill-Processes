%% Main script for running marker-in-cell method to simulate advection-dispersion processes

clear
% close all

%% Model parameters
%Initialize Solver (discretization of space and time)

% Spatial Discretization

% Define internodal boundaries (include at least all soillayer boundaries)
zBottom = -1;
dz = 0.05;
ModelDim = InitializeNodes('z', zBottom, dz);

alln = 1:ModelDim.znn;
allin = 1:ModelDim.znin;

% Boundary parameters
BoundaryPar.cTop = 0;
BoundaryPar.qTop = @qBoundary;
BoundaryPar.kSurf = 1e-2;
BoundaryPar.hAmb = -0.95;

% Time discretization
tRange = 0:0.2:15;

% Hydraulic properties
inFlow = -1e-1 ;                             % m^3 water per m^2 soil per day
d = 1e-3;                                    % diffusion coefficient 
theta = 0.4;                                 % volumetric water content (m^3 water/m^3 soil)
SoilPar.v(1:ModelDim.znin, 1) = inFlow;      % convective flux
SoilPar.d(1:ModelDim.znin, 1) = d;           % [m^2/s]
SoilPar.thetaN(1:ModelDim.znn, 1) = 0.3;     % m^3/m^3
SoilPar.thetaIn(1:ModelDim.znin, 1) = 0.3;   % m^3/m^3

hcap = 0.5; % Air entry head
SoilPar.alpha(alln, 1) = 1 ./ hcap;
SoilPar.n(alln, 1) = 1.9;
SoilPar.thetaR(alln, 1) = 0.04;
SoilPar.thetaS(alln, 1) = 0.4;
SoilPar.alphaIn(allin, 1) = 1 ./ hcap;
SoilPar.nIn(allin, 1) = 1.9;
SoilPar.thetaRIn(allin,1) = 0.04;
SoilPar.thetaSIn(allin,1) = 0.4;
SoilPar.kSatIn(allin,1) = 1e-2;

% Model initial states;
cIni(1:ModelDim.znn, 1)=1;

% Initialization of Markers and Cell approach:
% number of markers per cell
nMark = 20;
totMark = nMark * ModelDim.znn;

% Distribute markers over grid with a small random perturbation...
dzMark = (ModelDim.zin(end) - ModelDim.zin(1)) ./ totMark;
% rngSeed = 1;
% rng(rngSeed);
zMark(1:totMark, 1) = ((1:totMark) - 1/2) .* dzMark + dzMark .* (rand - 0.5); % the previous code allowed marker te be at x larger than zero, resulting in NaN's later in the script

% Define initial Concentration distribution for markers
cMark = interp1(ModelDim.zn, cIni, zMark, 'cubic', 'extrap');

%% RK4 approach with Matrix...
tic
difftime = 1e-11; % this is assumed to be zero (equality for time)
convcrit = 1e-11; % test value for assessing convergence

tMic = tRange(1);
tend = tRange(end);

% Max timestep is based on cell size and flow velocity
dtMax = max(abs(ModelDim.dzn)) ./ max(abs(SoilPar.v));

% Initialize output matrices
% store initial values in output matricese
timeOut(1) = 0;
coiTime(:, 1) = cIni;
c = cIni;
iTime = 1;

% This is to check mass balance
initialMass = ComputeSoluteMass(zMark, cMark, ModelDim);
massOutCumulative = 0;

% tic
% [qOut, thetaOut, hOut, time] = Richards(tRange, dtMax, ModelDim, SoilPar, BoundaryPar);
% inFlow = mean(mean(qOut));
% toc

tic

%% Simulation loop
while abs(tMic-tend) > difftime,
    fluidVelocity = SoilPar.v;
%     fluidVelocity = interp2(ModelDim.zin, time, qOut, ModelDim.zin, tMic, 'linear');
    [tMic, cMark, cNode, zMark, nMarkOut, massOut] = MarkerInCell(tMic, iTime, zMark, cMark, fluidVelocity, tRange, SoilPar, ModelDim, BoundaryPar);
    massOutCumulative = massOutCumulative + massOut;
    massRemaining = ComputeSoluteMass(zMark, cMark, ModelDim);
    % Update output matrices
    if abs(tMic - tRange(iTime + 1)) < difftime,
        iTime = iTime + 1;
        timeOut(iTime) = tMic;
        coiTime(:, iTime) = cNode;
    end
    
end
toc

tic
%% Running analytical solution
cAn = SoluteTransportAnalytic(ModelDim.zn, tRange, inFlow, d);
toc

%% Checking mass balance
finalMass = ComputeSoluteMass(zMark, cMark, ModelDim);
fprintf('Initially was: %f, flushed %f, remaining %f\n', initialMass, massOutCumulative, finalMass)
balance = abs(initialMass - massOutCumulative - finalMass);
fprintf('Balance is: %f\n', balance);

%% Plots

%Plot times series as a function of time for selected depths
figure(1)
clf
firstNodeIdx = 1;
displayStep = 5;
plot(timeOut, cAn(firstNodeIdx:displayStep:end, :), 'LineWidth', 2);
hold on
plot(timeOut, coiTime(firstNodeIdx:displayStep:end, :), 'x-', 'LineWidth', 1);
hold off
axBounds = axis;
axis([axBounds(1:2), 0, 1]);
xlabel('time [days]');
ylabel('Concentration [kg/m^3]');

clear