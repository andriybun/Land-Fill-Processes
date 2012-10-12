function y = Yvec(t,h,SoilPar,BoundaryPar,ModelDim)
    
    nn = ModelDim.nn;
    alln = 1:nn;
    nin = ModelDim.nin;
    allin = 1:nin;
    
    dzn = ModelDim.dzn;
    dzin = ModelDim.dzin;
    
    %Seff = SeffTHe(h,SoilPar);
    K =  HyC_function(h,SoilPar,ModelDim);
    %dKdh = dKdh_function(h,SoilPar,ModelDim);
    Ksrf = BoundaryPar.Ksurf;

    y(1,1) = +K(2,1)/dzin(1)+BoundaryPar.qTop(t)./dzin(1);
    
    y(2:nn-1,1) = (-K(2:nn-1,1) + K(3:nn,1))./dzin(2:nn-1,1);
 
    y(nn,1) = -Ksrf./dzin(nn).*BoundaryPar.hamb - K(nn)./dzin(nn);
end



