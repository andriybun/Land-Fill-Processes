function u = dKdh_function(h,SoilPar,ModelDim)
    
    nn = ModelDim.nn;
    nin = ModelDim.nin;
    dzin = ModelDim.dzin;
    
    alpha = SoilPar.alpha_in;
    n = SoilPar.n_in;
    m = 1-1./n;
    Ksat = SoilPar.Ksat_in;
    
    %Calculate internodal pressures
    hin(1,1) = h(1,1);
    hin(2:nn,1) = (dzin(1:nn-1,1).*h(1:nn-1,1) + dzin(2:nn,1).*h(2:nn,1))./...
        (dzin(1:nn-1,1)+dzin(2:nn));
    hin(nn+1,1) = h(nn,1);
    
    allin = 1:nin;
    
    Seff_in (allin,1) = (1+(alpha.*abs(hin)).^n).^-m.*(hin<0)+(hin>=0);
    
    %Note we only require dSeff/dh for the dK/dh, Ch includes the
    %compressibility of water, we need to subtract this term...
    %dSeffdh(allin,1) = Cm_function(hin,SoilPar)-Seff(allin1,1).*Sw;
    
    dSeffdh = m.*n.*alpha.*(alpha.*abs(hin)).^(n-1).*(1+(alpha.*abs(hin)).^n).^(-m-1) .* ...
        (hin<0)+0.*(hin>=0);
    
    %This model assumes a simplified krel = Seff^3
    u = 3.*Seff_in.^2.*dSeffdh.*Ksat;
end