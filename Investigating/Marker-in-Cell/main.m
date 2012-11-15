%% Template for a Simulation Model
%
% Example template to illustrate an algorithmic approach to running simulations
% in Matlab

clear
close all

%% Discretization
%Initialize Solver (discretization of space and time)

% Spatial Discretization
% x,y,z
% Define internodal boundaries (include at least all soillayer boundaries)
zin = [0:-0.05:-1]';
nin = length(zin);
zn(1,1) = zin(1,1);
zn(2:nin-2,1) = (zin(2:nin-2)+zin(3:nin-1))./2;
zn(nin-1,1) = zin(nin);
nn = length(zn);

dzn = zn(2:end,1)-zn(1:end-1,1);
dzin = zin(2:end,1)-zin(1:end-1,1);

ModelDim.zn=zn;
ModelDim.zin=zin;
ModelDim.dzn=dzn;
ModelDim.dzin=dzin;
ModelDim.nn=nn;
ModelDim.nin=nin;

% boundary parameters
BoundaryPar.CTop = 0;                   % Average yearly Concentration (kg/m³)

% Time Discretization
% start time
trange = 0:1:80;                        % chaged from 5 years
%% Model initialization

% model parameters
InFlow = -0.01 ;                        % m3 water per m2 soil per day
theta = 0.4;                            % volumetric water content (m3 water/m3 soil)
SoilPar.v(1:ModelDim.nin,1)= InFlow;    % convective flux
SoilPar.D(1:ModelDim.nin,1)= 0.0001;    % [m^2/s]
SoilPar.thetan(1:ModelDim.nn,1) = 0.3;  % m^3/m^3
SoilPar.thetain(1:ModelDim.nin,1)=0.3;  % m^3/m^3

% model initial states;
Cini(1:nn,1)=1;

%Initialization of Markers and Cell approach:
% number of markers per cell
nMark = 20;
totMark = nMark * nin;

%distribute markers over grid with a small random perturbation...
dzmrk = (zin(end) - zin(1)) ./ totMark;
zmrk(1:totMark, 1) = ((1:totMark) - 1/2) .* dzmrk + dzmrk .* (rand - 0.5); % the previous code allowed marker te be at x larger than zero, resulting in NaN's later in the script

%Define initial Concentration distribution for markers
Cmrk = interp1(zn, Cini, zmrk, 'cubic', 'extrap');

%% RK4 approach with Matrix...
tic
difftime = 1e-11; % this is assumed to be zero (equality for time)
convcrit = 1e-11; % test value for assessing convergence

t = trange(1);
tend = trange(end);

%max timestep is based on cell size and flow velocity
dtmax = 0.2.*max(abs(dzn))./max(abs(SoilPar.v));

%Initialize output matrices
%store initial values in output matricese
Time2(1) = 0;
Con2(1,:) = Cini';
C = Cini;
nout = 1;

%% Diff notation
cMark = Cmrk;
zMark = zmrk;
tMac = t;
%%

while abs(t-tend) > difftime,
    %%
    tMac = t;
    [tMac, cMark, cNode, zMark, nMarkOut] = mic(tMac, nout, zMark, cMark,  SoilPar.v, trange, SoilPar, ModelDim, BoundaryPar);
    C = cNode;
    t = tMac;
    %%
    
    %Update output matrix
    if abs(t-trange(nout+1)) < difftime,
        nout = nout+1;
        Time2(nout) = t;
        Con2(nout,:) = C';
    end
    
end
%Calculate Top Boundary condition for output
 %CBnd = BoundaryPar.CBottom;
 %Con2(:,1) = Con2';
toc



%% Close simulation

%Plot times series as a function of time for selected depths
figure(1)
clf
hold on
plot(Time2,Con2(:,:),'o-')

xlabel ('time [days]')
ylabel ('Concentration [kg/m^3]')

%plot temperature as a function of depth for selected times
figure(2)
clf
hold on
plot(Con2(:,:),zn,'o-')

xlabel ('Concentration [kg/m^3]')
ylabel ('depth [m]')
%save 'Solute_Advection_Diffusion01 '