function [ORI Meas Comp S Rp I] = bioreactorP()

% load stoichiometry, initial conditions and rate parameters from Excell
a = fopen('matrices_MSWS2.csv');
N = 13; M = 24; nr = 3+1; nc = 4; nout = 14; I.nCx = 3; XX = [];

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
Rp.Vlini = str2double(XX(4+3*nr+8,3));
Rp.rVl = str2double(XX(4+3*nr+8,4));

% Definition of Gas parameters
R = 0.082; % L*atm/K/mol 
T = 298.15; % K
p = 1; % atm
H_CO2 = 29.41; % L*atm/mol 
Vg = str2double(XX(4+3*nr+8,5)); % L
I.Vg = Vg;
Rp.rN2in = str2double(XX(4+3*nr+8,2)); % mol/L

% The initial amount of gasses are dependent on the initial pH and pCO2
Ci = initialize_MSWS2(Comp,1);  
k1 = find(strcmp('H2CO3.con',Comp.all));
I.CO2i = Ci(k1)*H_CO2*Vg/R/T; % mol px = Hx*Cx nx = pxVg/RT
k1 = find(strcmp('NH3.con',Comp.all));
I.NH3i = Ci(k1)*0.0164*Vg/R/T;
I.H2Oi = 0.0316*Vg/R/T;
I.N2i = p*Vg/R/T-I.CO2i-I.NH3i-I.H2Oi;
I.CH4i = 0;
I.H2Si = 0;

% The initial amount of H2CO3.tot is dependent on the initial pCO2
k1 = find(strcmp('H2CO3.tot',Comp.all));
Comp.masteri(k1) = Ci(k1);  % set initial H2CO3.tot

fclose('all');

% initialize ORCHESTRA
ORI = initialize_MSWS2(Comp,0); 
   
% Import experimental data 
a = fopen('Data_MSWS2.csv');
N = 14; M = 225; XX2 = [];

for i=1:M
    X = textscan(a,'%s',N);
    XX2 = [XX2;X{:}'];
end

Biom = str2double(XX2(2:225,2))/24.5/1000;
Biom(216) = NaN;
pHm = str2double(XX2(2:225,5));
pHm(45) = NaN;
VFAm = str2double(XX2(2:58,7))/1000/60.05;
tm = str2double(XX2(2:225,1));  % t --> Biom and pHm
tm2 = str2double(XX2(2:58,6)); % t --> VFAm
tm3 = str2double(XX2(2:132,8)); % t --> NH4m
NH4m = str2double(XX2(2:132,9))/1000/18;  
ECm = str2double(XX2(2:225,12)); % mS/cm
CH4pm = str2double(XX2(2:225,13))/100; % percentage
CO2pm = str2double(XX2(2:225,14))/100; % percentage

fclose('all');

Meas.tm = [tm;tm;tm2];  
Meas.data = [Biom;pHm;VFAm];
kNaN = find(isnan(Meas.data));
Meas.tm(kNaN) = [];
Meas.data(kNaN) = [];
Meas.ldata = [length(Biom) length(pHm) length(VFAm)];

end