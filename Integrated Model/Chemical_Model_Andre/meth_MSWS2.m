function dCTdt = meth_MSWS2(t, CT, S, Comp, Rp, ORI, rVl, nCx, Vg, rN2in, F)

[nreac,ncomp] = size(S);
V = CT(end-11);
N2 = CT(end-10);
CO2 = CT(end-9);
CH4 = CT(end-8);
H2S = CT(end-7);
NH3 = CT(end-6);
H2O = CT(end-5);
Ca = CT(end-2);
Na = CT(end-1);
Cl = CT(end);
CT = CT(1:end-12);

% kick out negative concentrations to stabilize ORCHESTRA calculations
for i=1:ncomp
    if CT(i) < 1e-9 || CT(i) > 100000
        if i ~= find(strcmp('H+.tot',Comp.master))
           CT(i) = 1e-9;
        end
    end
end

% calculate speciation concentrations from total concentrations
global C C2 Vv
n = (1:1:(length(CT)+length(Comp.consti)));
CT1 = [CT;Comp.consti']/V;
ro = ORI.Calculate(n, CT1); 
k = length(Comp.extrai);
C = ro(length(CT)+length(Comp.consti)+1:end-k);
C2 = ro(end-length(Comp.extrai)+1:end);
Vv = V;

% Calculate substrate limitation factor for each master species
Ks = Rp.Ks;
fs = ones(nreac,1);
for ii = 1:nreac,
  idx = find(Ks(ii,:)>0);
  fs(ii,1) = prod(CT(idx)'./(CT(idx)' + Ks(ii,idx))); % total acid concentration serves as substrate
%   fs(ii,1) = prod(C(idx)'./(C(idx)' + Ks(ii,idx))); % in case that only protonated acids are substrate
end

% Calculate pH limitation factor for each reaction
pH = C(find(strcmp('H+.tot',Comp.master)));
fpH = ones(nreac,1); 
k1 = find(Rp.pH1 ~= 0);
pH1 = Rp.pH1(k1);
pH2 = Rp.pH2(k1);
fpH(k1) = smoothstep(pH1,pH2,pH);

% Calculate Temperature/Ammonia limitation factor for each reaction
fT = ones(nreac,1);
k1 = find(Rp.NH3i ~= 0);
fT(k1) = 1/(1+(C2(find(strcmp('NH3.con',Comp.extra)))/Rp.NH3i(k1)));

% Calculate total limitation factor for each reaction
ftot = fs.*fpH.*fT;

% Calculate uptake rate for each each reaction (mol-S/C-mol-X/d)
Rmax = Rp.max{2};
R = Rmax.*ftot;

% Correction of uptake rate for reactions that compete for limiting substrate
k1 = find(strcmp('hydrolysis',Rp.max{1}));
k2 = find(strcmp('Desulfurication',Rp.max{1}));
ff = 1;
Th3 = R(k1)/(R(k1)+ff*R(k2));
Th6 = 1-Th3;
R(k1) = R(k1)*Th3;
R(k2) = R(k2)*Th6;

% uptake rate matrix  (uptake rate per reaction per master species)
B0 = (S.*(R*ones(1,ncomp)));   

% process matrix (ncompxncomp --> uptake rate per master species per Cx or Cs)
Bpx = zeros(ncomp,ncomp); 
Cxb = ncomp+1-nCx;
Cxe = ncomp;
for i = 1:ncomp
    for j = 1:nreac
        if S(j,1)==-1;
            Bpx(i,1) = Bpx(i,1) + B0(j,i);
        else
            kk = find(S(j,Cxb:Cxe)~=0)+Cxb-1;
            if kk ~= 0
                Bpx(i,kk) = Bpx(i,kk) + B0(j,i);
            else
                kk = find(S(j,:)==-1);
                Bpx(i,kk) = Bpx(i,kk) + B0(j,i);
            end
        end
    end
end

Rt = Bpx*CT;   % mol/day

% calculation diffusion transport term
F(1:length(Comp.master),1)=F;
F(5:6) = 0;
Rt = Rt+F.*CT;
F(7) = 2*Ca*F(7)+(Na-Cl)*F(7);
Rt(7) = Rt(7)-F(7);

% calculation of rates of gas phases
k1 = find(strcmp('H2CO3.con',Comp.out));
k2 = find(strcmp('H2CO3.tot',Comp.master));
k3 = find(strcmp('CH4_(cum)',Comp.master));
k4 = find(strcmp('H2S.tot',Comp.master));
k5 = find(strcmp('H2S.con',Comp.extra));
k6 = find(strcmp('NH3.tot',Comp.master));
k7 = find(strcmp('NH3.con',Comp.extra));
k8 = find(strcmp('H2O.tot',Comp.master));

Rin_CO2 = -100*(((CO2*0.082*298.15/Vg)/29.41)-C(k1))*V; % mol/d
Rt(k2) = Rt(k2)-Rin_CO2;
Rin_H2S = -100*(((H2S*0.082*298.15/Vg)/10)-C2(k5))*V; % mol/d
Rt(k4) = Rt(k4)-Rin_H2S;
Rin_NH3 = -100*(((NH3*0.082*298.15/Vg)/0.0164)-C2(k7))*V; % mol/d
Rt(k6) = Rt(k6)-Rin_NH3;
Rin_H2O = -100*(((H2O*0.082*298.15/Vg)/0.00056835)-55.6);
Rt(k8) = Rt(k8)-Rin_H2O;

% addition of volume and gas phases to dx/dt vector
% dV
Rt = [Rt;rVl]; % rV
% dN2
rN2 = rN2in-(rN2in+Rin_CO2+Rt(k3)+Rin_H2S+Rin_NH3+Rin_H2O)*N2/(N2+CO2+CH4+H2S+NH3+H2O);
Rt = [Rt;rN2]; % rCO2g
% dC02
rCO2 = Rin_CO2-(rN2in+Rin_CO2+Rt(k3)+Rin_H2S+Rin_NH3+Rin_H2O)*CO2/(N2+CO2+CH4+H2S+NH3+H2O);
Rt = [Rt;rCO2]; % rCO2g
% dCH4
rCH4 = Rt(k3)-(rN2in+Rin_CO2+Rt(k3)+Rin_H2S+Rin_NH3+Rin_H2O)*CH4/(N2+CO2+CH4+H2S+NH3+H2O);
Rt = [Rt;rCH4]; % rCH4g
% dH2S
rH2S = Rin_H2S-(rN2in+Rin_CO2+Rt(k3)+Rin_H2S+Rin_NH3+Rin_H2O)*H2S/(N2+CO2+CH4+H2S+NH3+H2O);
Rt = [Rt;rH2S]; % rH2Sg
% dNH3
rNH3 = Rin_NH3-(rN2in+Rin_CO2+Rt(k3)+Rin_H2S+Rin_NH3+Rin_H2O)*NH3/(N2+CO2+CH4+H2S+NH3+H2O);
Rt = [Rt;rNH3]; % rNH3g
% dH2O
rH2O = Rin_H2O-(rN2in+Rin_CO2+Rt(k3)+Rin_H2S+Rin_NH3+Rin_H2O)*H2O/(N2+CO2+CH4+H2S+NH3+H2O);
Rt = [Rt;rH2O]; % rH2Og

% addition of cumulative CO2 and CH4 to dx/dt vector
rCO2_out = (rN2in+Rin_CO2+Rt(k3)+Rin_H2S+Rin_NH3+Rin_H2O)*CO2/(N2+CO2+CH4+H2S+NH3+H2O);
Rt = [Rt;rCO2_out]; 
rCH4_out = (rN2in+Rin_CO2+Rt(k3)+Rin_H2S+Rin_NH3+Rin_H2O)*CH4/(N2+CO2+CH4+H2S+NH3+H2O);
Rt = [Rt;rCH4_out];

% addition of mass transfer of constant components
Rt = [Rt;F(1)*Ca;F(1)*Na;F(1)*Cl];

dCTdt = Rt; % R*(Cx or Cs) per reaction (mol-S/L/d)
end
