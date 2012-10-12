function S = Kmat(h,SoilPar,BoundaryPar,ModelDim)
    
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

%     wf1 = 0;
%     wf2 = 0;
%     wf3 = dzin(1)./(dzin(1)+dzin(2));
%     wf4 = dzin(2)./(dzin(1)+dzin(2));
%     
%     %upper boundary
%     a(1,1) = 0;
%     b(1,1) = -(K(2,1) + h(1,1).*dKdh(2,1).*wf3)./ ...
%         (dzin(1).*dzn(1))+dKdh(2,1).*wf3./dzin(1);
%     c(1,1) = (K(2,1)+h(2,1).*dKdh(2,1).*wf4)./ ...
%         (dzin(1,1).*dzn(1,1)) + dKdh(2,1).*wf4./dzin(1);
%      
%     %middel nodes
%      
%     wf1 = dzin(1:nn-2)./(dzin(1:nn-2)+dzin(2:nn-1));
%     wf2 = dzin(2:nn-1)./(dzin(1:nn-2)+dzin(2:nn-1));
%     wf3 = dzin(2:nn-1)./(dzin(2:nn-1)+dzin(3:nn));
%     wf4 = dzin(3:nn)./(dzin(2:nn-1)+dzin(3:nn));
%     
%     a(2:nn-1,1) = (K(2:nn-1,1) + h(1:nn-2,1).*dKdh(2:nn-1,1).*wf1)./ ...
%         (dzin(2:nn-1).*dzn(1:nn-2)) - dKdh(2:nn-1,1).*wf1./dzin(2:nn-1);
%     
%     b(2:nn-1,1) = -(K(2:nn-1,1)+h(2:nn-1,1).*dKdh(2:nn-1,1).*wf2)./ ...
%         (dzin(2:nn-1,1).*dzn(1:nn-2)) - ...
%         (K(3:nn,1)+h(2:nn-1,1).*dKdh(3:nn,1).*wf3)./ ...
%         (dzin(2:nn-1,1).*dzn(2:nn-1,1)) - ...
%         (dKdh(2:nn-1,1).*wf2-dKdh(3:nn,1).*wf3)./dzin(2:nn-1,1);
%      
%     c(2:nn-1,1) = (K(3:nn,1)+h(3:nn,1).*dKdh(3:nn,1).*wf4)./ ...
%         (dzin(2:nn-1).*dzn(2:nn-1,1)) + dKdh(3:nn,1).*wf4./dzin(2:nn-1);
%      
%     %lower boundary
%     wf1 = dzin(nn-1)./(dzin(nn-1)+dzin(nn));
%     wf2 = dzin(nn)./(dzin(nn-1)+dzin(nn));
%     wf3 = 0;
%     wf4 = 0;
%     
%     a(nn) = (K(nn)+h(nn-1,1).*dKdh(nn).*wf1)./ ...
%         (dzin(nn).*dzn(nn-1)) - dKdh(nn).*wf1./dzin(nn);
%     b(nn) = -(K(nn)+h(nn,1).*dKdh(nn).*wf2)./ ...
%         (dzin(nn).*dzn(nn-1)) + ...
%         Ksrf./dzin(nn) - dKdh(nn).*wf2./dzin(nn);
%     c(nn) = 0;
     
    a(1,1) = 0;
    b(1,1) = -K(2,1)./ ...
        (dzin(1,1).*dzn(1,1));
    c(1,1) = K(2,1)./ ...
        (dzin(1,1).*dzn(1,1));
    
    %middel nodes
    a(2:nn-1,1) = K(2:nn-1,1)./ ...
        (dzin(2:nn-1).*dzn(1:nn-2));
    
    b(2:nn-1,1) = -K(2:nn-1,1)./ ...
        (dzin(2:nn-1,1).*dzn(1:nn-2)) - ...
        K(3:nn,1)./ ...
        (dzin(2:nn-1,1).*dzn(2:nn-1,1));
    
    c(2:nn-1,1) = K(3:nn,1)./ ...
        (dzin(2:nn-1).*dzn(2:nn-1,1));
    
    %lower boundary
    a(nn) = K(nn)./ ...
        (dzin(nn).*dzn(nn-1));
    b(nn) = -K(nn)./ ...
        (dzin(nn).*dzn(nn-1)) + ...
        Ksrf./dzin(nn);
    c(nn) = 0;
    
    
    
    B =  diag(a(2:nn),-1)+diag(b,0)+diag(c(1:nn-1),+1);
    S = sparse(B);
end



