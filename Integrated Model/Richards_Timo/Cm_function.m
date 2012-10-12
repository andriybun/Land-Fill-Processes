function Cm=Cm_function(h,SoilPar)
alpha = SoilPar.alpha;% m^-1 Air entry value
n = SoilPar.n;%van Genuchten parameter
m = 1-1./SoilPar.n;%van Genuchten parameter
%Ksat = SoilPar.Ksat;%m/s saturated hydraulic conductivity
thetaS = SoilPar.thetaS;%Saturated water content
thetaR = SoilPar.thetaR;%Residual water content

rhow=1000;
g = 9.81;%m/s^2 accleration due to gravity
Sw = 4e-10.*rhow.*g; % compressibility of water

Seff = (1+(alpha.*abs(h)).^n).^-m.*(h<0)+(h>=0);
S = Seff+(thetaR./thetaS);
C = m.*n.*alpha.*(alpha.*abs(h)).^(n-1).*(1+(alpha.*abs(h)).^n).^(-m-1) .* ...
   (h<0)+0.*(h>=0);
C = C + Sw.*S;
Cm = sparse(diag(C,0));
end

