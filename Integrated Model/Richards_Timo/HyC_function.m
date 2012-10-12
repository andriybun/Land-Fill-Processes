function u = HyC_function(h,SoilPar,ModelDim)
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
    krw = Seff_in.^3;
    u = Ksat .* krw;
      
end

