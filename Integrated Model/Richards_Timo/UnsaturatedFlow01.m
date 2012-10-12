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
BoundaryPar.qTop = -0.002; %.
BoundaryPar.Ksurf = 1e-2;
BoundaryPar.hamb = -0.1;

% Time Discretization
% start time
trange = [0 250];%0:5:250*365; %5 years)
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
hbot = -0.95;
hIni(alln,1)= hbot-(zn-zref);

%% Dynamic loop (non-stationary model) %% Implicit solution with ode15i
%options = odeset('RelTol',1e-5,'Stats','on',...
%   'Mass',@(t,h) Cm_function(h,SoilPar),'MStateDependence','strong');


options = odeset('AbsTol',1e-5,'RelTol',1e-5,'Stats','on',...
    'Jacobian', @(t,h) JacobianRichards(h,SoilPar,BoundaryPar,ModelDim),...
    'Mass',@(t,h) Cm_function(h,SoilPar),'MStateDependence','strong');



tic

[Time,hSim] = ode23t(@(t,h) NettoFlux(t,h,SoilPar,ModelDim,BoundaryPar),...
    trange,hIni,options);

toc

figure(1)
clf
hold on
plot(hSim,zn);
%plot temperature as a function of depth for selected times
figure(2)
clf
hold on
plot(Time,hSim)

thetaOut = [];
qOut = [];
%%Nested functions for Jacobian and Massmatrix
for ii=1:length(Time),
    [States] = Waterstates(Time(ii),hSim(ii,:)',SoilPar,ModelDim,BoundaryPar);
    thetaOut = [thetaOut;States.theta'];
    qOut = [qOut;States.q'];
end
MassBalance
%end
