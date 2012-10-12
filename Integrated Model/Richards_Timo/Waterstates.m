function [States] = Waterstates(Time,h,SoilPar,ModelDim,BoundaryPar)
    
    Ksurf=BoundaryPar.Ksurf;
    hamb=BoundaryPar.hamb;
    %Ksat = repmat(SoilPar.Ksat_in(:)',nt,1);
    
    alpha = SoilPar.alpha;
    n = SoilPar.n;
    m = 1-1./n;
    thetaR = SoilPar.thetaR;
    thetaS = SoilPar.thetaS;
    
    Seff = (1+(alpha.*abs(h)).^n).^-m.*(h<0)+(h>=0);
    theta = thetaR + Seff.*(thetaS-thetaR);
    
    q = Richards1DTHe(Time,h,SoilPar,ModelDim,BoundaryPar);
    States.Seff = Seff;
    States.theta=theta;
    States.q=q;
