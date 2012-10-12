%Calculate MassBalance of Simulation
%Data: States.theta water content (ntime, nnodes)
%      States.q     waterflux (ntime, nnodes+1)

%Total change of mass in column = spatially integrated water content the final timestep -
%spatially integrated water content at first timestep
dt = Time(2:end)-Time(1:end-1);

%unit: m3/m2;
TotWatMass = -thetaOut(:,1:nn)*dzin; 
NFlux = -(qOut(:,1)-qOut(:,end));
% Total change of mass in column should be equal to total the difference 
% between the time integrated flux at the top and the bottom boundaries

MassBal1 = diff(TotWatMass)-NFlux(2:end).*dt;
MassBal2 = diff(TotWatMass)-NFlux(1:end-1).*dt;
MassBal3 = gradient(TotWatMass,Time)-NFlux;


figure(3)
plot(Time(2:end),MassBal1,'b');
hold on
plot(Time(2:end),MassBal2,'g');
%plot(Time,MassBal3,'m')
