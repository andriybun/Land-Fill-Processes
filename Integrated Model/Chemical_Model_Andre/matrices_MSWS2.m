%-----------------------------------------------------Implementation of Reichel+Kleerebezem+Meeussen-------------------------------------------    

%-----------------------------------------------------------------Defining-matrices------------------------------------------------------------
clc
clear all
close all

pname = './';
addpath(genpath(pname));

% load stoichiometry, initial conditions and rate parameters from Excell
a = fopen('matrices_MSWS2.csv');
N = 13; M = 24; nr = 3+1; nc = 4; nout = 14; nCx = 3; XX = [];

for i=1:M
    X = textscan(a,'%s',N);
    XX = [XX;X{:}'];
end

% Definition of Chemical equilibrium system + initial conditions ((C-)mol/L)
Comp.master = XX(1,2:end);
Comp.masteri = str2double(XX(2,2:end));

Comp.out = XX(3,2:end);
Comp.outi = str2double(XX(4,2:end)); 

Comp.const = XX(4+3*nr+1,2:nc+1);
Comp.consti = str2double(XX(4+3*nr+2,2:1+nc));

Comp.extra = [XX(4+3*nr+3,2:nout/2+1) XX(4+3*nr+5,2:nout/2+1)];
Comp.extrai = str2double([XX(4+3*nr+4,2:nout/2+1) XX(4+3*nr+6,2:nout/2+1)]);

Comp.all = [Comp.master Comp.const Comp.out Comp.extra];
Comp.alli = [Comp.masteri Comp.consti Comp.outi Comp.extrai];

% Definition stoichiomatrix
S = str2double(XX(4+2:4+nr,2:end));

% Definition of half saturation constants (mol/L) 
Rp.Ks = str2double(XX(4+2+nr:4+2*nr,2:end));

% Definition maximum uptake rates (mol-S/C-mol-X/d)
Rp.max = {XX(4+2+2*nr:4+3*nr,1) str2double(XX(4+2+2*nr:4+3*nr,2))};

% Definition parameters pH inhibition
Rp.pH1 = str2double(XX(4+2+2*nr:4+3*nr,3));
Rp.pH2 = str2double(XX(4+2+2*nr:4+3*nr,4));

% Definition parameters ammonia inhibition
Rp.NH3i = str2double(XX(4+2+2*nr:4+3*nr,5));

% Definition of Volume parameters
Vlini = str2double(XX(4+3*nr+8,3));
rVl = str2double(XX(4+3*nr+8,4));

% Definition of Gas parameters
R = 0.082; % L*atm/K/mol 
T = 298.15; % K
p = 1; % atm
H_CO2 = 29.41; % L*atm/mol 
Vg = str2double(XX(4+3*nr+8,5)); % L
rN2in = str2double(XX(4+3*nr+8,2)); % mol/L

% The initial amount of gasses are dependent on the initial pH and pCO2
Ci = initialize_MSWS2(Comp,1);  
k1 = find(strcmp('H2CO3.con',Comp.all));
CO2i = Ci(k1)*H_CO2*Vg/R/T; % mol px = Hx*Cx nx = pxVg/RT
k1 = find(strcmp('NH3.con',Comp.all));
NH3i = Ci(k1)*0.0164*Vg/R/T;
H2Oi = 0.0316*Vg/R/T;
N2i = p*Vg/R/T-CO2i-NH3i-H2Oi;
CH4i = 0;
H2Si = 0;

% The initial amount of H2CO3.tot is dependent on the initial pCO2
k1 = find(strcmp('H2CO3.tot',Comp.all));
Comp.masteri(k1) = Ci(k1);  % set initial H2CO3.tot

% Definition of transport parameter
F = str2double(XX(4+3*nr+8,7)); % diffusion parameter 1/day

fclose('all');
%----------------------------------------------------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------Simulation---------------------------------------------------------------
% initialize ORCHESTRA
ORI = initialize_MSWS2(Comp,0); 

% global variables to store intermediate ORCHESTRA results
global Call Call2 V2 tt
Call = [];  Call2 = []; V2 = []; tt = []; 

% Solve Differential Equation
CO2_outi = 0;
CH4_outi = 0;
Const = Comp.consti(2:4);
CTini = [Comp.masteri*Vlini Vlini N2i CO2i CH4i H2Si NH3i H2Oi CO2_outi CH4_outi Const];
Comp.consti = Comp.consti*Vlini;
options = odeset('OutputFcn','StoreC','Refine',1, 'AbsTol', 1e-12, 'RelTol', 1e-7);
[t,CT] = ode15s(@meth_MSWS2, [0:1:740], CTini, options, S, Comp, Rp, ORI, rVl, nCx, Vg, rN2in, F);

%----------------------------------------------------------------------------------------------------------------------------------------------

%-------------------------------------------------------------------Data-processing------------------------------------------------------------

% load experimental data 
a = fopen('Data_MSWS2.csv');
N = 14; M = 225; XX2 = [];

for i=1:M
    X = textscan(a,'%s',N);
    XX2 = [XX2;X{:}'];
end

Biom = str2double(XX2(2:225,2))/24.5/1000;
pHm = str2double(XX2(2:225,5));
VFAm = str2double(XX2(2:58,7))/1000/60.05;
tm = str2double(XX2(2:225,1));  % t --> Biom and pHm
tm2 = str2double(XX2(2:58,6)); % t --> VFAm
tm3 = str2double(XX2(2:132,8)); % t --> NH4m
NH4m = str2double(XX2(2:132,9))/1000/18;  
ECm = str2double(XX2(2:225,12)); % mS/cm
CH4pm = str2double(XX2(2:225,13))/100; % percentage
CO2pm = str2double(XX2(2:225,14))/100; % percentage

fclose('all');

% Define simulated data
V = CT(:,end-11);
N2 = CT(:,end-10);
CO2 = CT(:,end-9);
CH4 = CT(:,end-8);
H2S = CT(:,end-7);
NH3 = CT(:,end-6);
H2O = CT(:,end-5);
CO2_out = CT(:,end-4);
CH4_out = CT(:,end-3);
CT = CT(:,1:end-12);
CTl = length(CT(1,:));
ng = nCx;
Bios = (CO2_out+CH4_out)/1000;
k1 = find(strcmp('H[Acetate].tot',Comp.master));
VFAs = CT(:,k1)./V;
k1 = find(strcmp('pH',Comp.out));
pHs = Call(k1,:);

% Display development of Master Species in time (Total concentrations)
figure(1)
CT(:,6) = CT(:,6)/55.6*Vlini; % L water 

for ii = 1:CTl-ng,
  subplot(3,4,ii)
  hold on
  plot(t,CT(:,ii)/Vlini,'-r');
  title (Comp.master{ii});
end

subplot(3,4,CTl-2)
hold on
plot(t,CT(:,CTl-2),'-r');
plot(t,CT(:,CTl-1),'-b');
plot(t,CT(:,CTl),'-g');
title ('Biomass');
legend(Comp.master(CTl-2:CTl),'Location','BestOutside');

% Display development of concentrations in time 
figure(2)
for ii = 1:length(Call(:,1)),
  subplot(4,4,ii)
  hold on
  plot(tt,Call(ii,:),'-b');
  title (Comp.out{ii});
  %pause (0.01);
end

subplot(4,4,6)
plot(tt,Call2(14,:),'-b');
title('Calcite');

subplot(4,4,13)
plot(t,CO2_out,'-b');
title('CO2(cum)');
axis([0 730 0 1200])

subplot(4,4,14)
plot(t,CH4_out,'-b');
title('CH4(cum)');
axis([0 730 0 1200])

subplot(4,4,15)
plot(tt,Call2(2,:),'-b');
title(Comp.extra{2});

subplot(4,4,16)
plot(tt,Call2(2,:)+Call(3,:),'-b');
title('H2CO3 + HCO3');

% Displays model results with experimental data (biogascum, VFA, pH)
figure(3)
subplot(3,1,1)
hold on
plot(t,Bios,'-xb')
plot(tm,Biom,'xr')
title('Cumulative CO2+CH4 (kmol)')
xlabel('Time (days)')
ylabel('CO2+CH4 (kmol)')

subplot(3,1,2)
plot(t,VFAs,'-xb')
hold on
plot(tm2,VFAm,'xr' )
title('Total VFA (mol/L)')
xlabel('Time (days)')
ylabel('VFA (mol/L)')

subplot(3,1,3)
plot(tt,pHs,'-xb')
hold on 
plot(tm,pHm,'xr')
title('pH')
xlabel('Time (days)')
ylabel('pH')

% Calculate EC
EC_Ca = 119*Call2(1,:);         EC_HCO3 = 44.5*Call2(2,:);     EC_Ac = 38.6*Call2(3,:);
EC_CaAc = 119/2/2*Call2(4,:);   EC_Cl = 76.3*Call2(5,:);       EC_CO3 = 44.5*2*Call2(6,:);
EC_NaCO3 = 44.5/2*Call2(8,:);   EC_OH = 198.6*Call2(9,:);      EC_NaSO4 = 80/2/2*Call2(10,:);
EC_HSO4 = 52*Call2(11,:);       EC_Na = 50.1*Call2(12,:);      EC_NH4 = 73.5*Call(4,:);
EC_H = 349.8*(10.^-pHs);        EC_SO4 = 80*Call(8,:)   ;      EC_HS = 65*Call(9,:);

EC = EC_Ca + EC_HCO3 + EC_Ac + EC_CaAc + EC_Cl + EC_CO3 + EC_NaCO3 + EC_OH + EC_NaSO4 + EC_HSO4 + EC_Na + EC_NH4 + EC_H + EC_SO4 + EC_HS; % mS/cm 
 
% Displays EC per species 
% figure(4)
% plot(tt,EC_Ca,'-b')
% hold on
% plot(tt,EC_HCO3,'-r')
% plot(tt,EC_Ac,'-y')
% plot(tt,EC_CaAc,'-k')
% plot(tt,EC_Cl,'-m')
% plot(tt,EC_CO3,'-g')
% plot(tt,EC_NaCO3,'-p')
% plot(tt,EC_OH,'xb')
% plot(tt,EC_NaSO4,'xr')
% plot(tt,EC_HSO4,'xy')
% plot(tt,EC_Na,'xk')
% plot(tt,EC_NH4,'xm')
% plot(tt,EC_H,'xg')
% plot(tt,EC_SO4,'ob')
% plot(tt,EC_HS,'ok')

% Display modelled and experimental NH4, EC, pCH4 and pCO2
figure(5)
subplot(2,2,1)
plot(tm3,NH4m,'xr')
hold on
plot(tt,Call(4,:),'-xb')
title('NH4')
xlabel('Time (days)')
ylabel('NH4')

subplot(2,2,2)
plot(tm,ECm,'xr')
hold on
plot(tt,EC,'-xb')
title('EC')
xlabel('Time (days)')
ylabel('EC')

subplot(2,2,3)
plot(tm,CH4pm,'xr')
hold on
plot(t,CH4*0.082*298.15/80,'-xb')
title('%CH4')
xlabel('Time (days)')
ylabel('%CO2')

subplot(2,2,4)
plot(tm,CO2pm,'xr')
hold on
plot(t,CO2*0.082*298.15/80,'-xb')
title('%CO2')
xlabel('Time (days)')
ylabel('%CO2')

% Electroneutrality check
Ca = 2*Call2(1,:).*V2;       HCO3 = -1*Call2(2,:).*V2   ;    Ac = -1*Call2(3,:).*V2   ;
CaAc = 1*Call2(4,:).*V2;     Cl = -1*Call2(5,:).*V2   ;      CO3 = -2*Call2(6,:).*V2   ;
NaCO3 = -1*Call2(8,:).*V2;   OH = -1*Call2(9,:).*V2   ;      NaSO4 = -1*Call2(10,:).*V2   ;
HSO4 = -1*Call2(11,:).*V2;   Na = 1*Call2(12,:).*V2   ;      NH4 = 1*Call(4,:).*V2   ;
H = 1*(10.^-pHs).*V2;        SO4 = -2*Call(8,:).*V2   ;      HS = -1*Call(9,:).*V2;          

EN = Ca+HCO3+Ac+CaAc+Cl+CO3+NaCO3+OH+NaSO4+HSO4+Na+NH4+H+SO4+HS;

%------------------------------------------------------------------------------------------------------------------------------------------------