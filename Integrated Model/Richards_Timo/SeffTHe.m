function Seff = SeffTHe(h,SoilPar)
alpha = SoilPar.alpha;% m^-1 Air entry value
n = SoilPar.n;%van Genuchten parameter
m = 1-1./SoilPar.n;%van Genuchten parameter
Ksat = SoilPar.Ksat_in;%m/s saturated hydraulic conductivity
Seff = (1+(alpha.*abs(h)).^n).^-m.*(h<0)+(h>=0);
%krel = Seff.^0.5 .* (1-(1-Seff.^(1./m)).^m).^2.*(h<0)+(h>=0);
%u=krel*Ksat;
end
