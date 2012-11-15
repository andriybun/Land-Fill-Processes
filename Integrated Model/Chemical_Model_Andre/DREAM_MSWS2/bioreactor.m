function [sim] = bioreactor(P,ExtraPar)

S = ExtraPar.S; Rp = ExtraPar.Rp; Comp = ExtraPar.Comp; ORI = ExtraPar.ORI;
tm = ExtraPar.Meas.tm; I = ExtraPar.I; ldata = ExtraPar.Meas.ldata;
Vlini = Rp.Vlini;   rVl = Rp.rVl;   rN2in = Rp.rN2in; nCx = I.nCx;  Vg = I.Vg;

% global variables to store intermediate ORCHESTRA results
global Call Call2 V2 tt
Call = [];  Call2 = []; V2 = []; tt = []; 

% replacing parameters to be fitted by DREAM
Rpm = Rp.max{2}; Rpm(1) = P(1); Rpm(2) = P(2); Rp.max{2} = Rpm; 

% Solve Differential Equation
CO2_outi = 0;
CH4_outi = 0;
CTini = [Comp.masteri*Vlini Vlini I.N2i I.CO2i I.CH4i I.H2Si I.NH3i I.H2Oi CO2_outi CH4_outi];
Comp.consti = Comp.consti*Vlini;
options = odeset('OutputFcn','StoreC','Refine',1, 'AbsTol', 1e-12, 'RelTol', 1e-7);
[t,CT] = ode15s(@meth_MSWS2, [0:1:740], CTini, options, S, Comp, Rp, ORI, rVl, nCx, Vg, rN2in);

% Define simulated data
V = CT(:,end-8);
N2 = CT(:,end-7);
CO2 = CT(:,end-6);
CH4 = CT(:,end-5);
H2S = CT(:,end-4);
NH3 = CT(:,end-3);
H2O = CT(:,end-2);
CO2_out = CT(:,end-1);
CH4_out = CT(:,end);
CT = CT(:,1:end-9);
CTl = length(CT(1,:));
ng = nCx;
Bios = (CO2_out+CH4_out)/1000;
k1 = find(strcmp('H[Acetate].tot',Comp.master));
VFAs = CT(:,k1)./V;
k1 = find(strcmp('pH',Comp.out));
pHs = Call(k1,:);

% select simulated data that correspond to experimental data (same t) 
tma  = tm(1:ldata(1));
Bioa = zeros(length(tma),1); 
for i = 1:length(tma)
    k1 = find(t==tma(i));
    Bioa(i) = Bios(k1);
end

tma  = tm(ldata(1)+1:ldata(1)+ldata(2));
pHa = zeros(length(tma),1); 
for i = 1:length(tma)
    k1 = find(tt==tma(i));
    if k1 > 0
        pHa(i) = pHs(k1);
    else
        k1 = find(tt<tma(i));
            if k1 > 0;
                k1 = k1(end);
                diff = 1;
                pHav = (pHs(k1)+pHs(k1+diff))/2;
                pHa(i) = pHav;
            else
                k1 = find(tt>tma(i));
                k1 = k1(1);
                diff = 1;
                pHav = (pHs(k1)+pHs(k1+diff))/2;
                pHa(i) = pHav;
            end
    end
end

tma  = tm(ldata(1)+ldata(2)+1:end);
VFAa = zeros(length(tma),1);
for i = 1:length(tma)
    k1 = find(t==tma(i));
    if k1 > 0
        VFAa(i) = VFAs(k1);
    end
    
    k1 = find(t<tma(i));
    if k1 > 0;
        k1 = k1(end);
        diff = 1;
        VFAav = (VFAs(k1)+VFAs(k1+diff))/2;
        VFAa(i) = VFAav;
    else
        k1 = find(t>tma(i));
        k1 = k1(1);
        diff = 1;
        VFAav = (VFAs(k1)+VFAs(k1+diff))/2;
        VFAa(i) = VFAav;
    end
end

% Displays model results with experimental data (biogascum, VFA, pH) 
% figure(3)
% subplot(3,1,1)
% hold on
% plot(t,Bios)
% Biom = ExtraPar.Meas.data(1:length(tm));
% plot(tm,Biom,'xr')
% title('Cumulative CO2+CH4 (kmol)')
% xlabel('Time (days)')
% ylabel('CO2+CH4 (kmol)')
% 
% subplot(3,1,2)
% plot(t,VFAs)
% hold on
% VFAm = ExtraPar.Meas.data(2*length(tm)+1:end);
% plot(tm2,VFAm,'xr' )
% title('Total VFA (mol/L)')
% xlabel('Time (days)')
% ylabel('VFA (mol/L)')
% 
% subplot(3,1,3)
% plot(tt,pHs)
% hold on
% pHm = ExtraPar.Meas.data(length(tm)+1:2*length(tm));
% plot(tm,pHm,'xr')
% title('pH')
% xlabel('Time (days)')
% ylabel('pH')

sim = [Bioa;pHa;VFAa];

end