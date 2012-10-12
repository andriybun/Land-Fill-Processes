function NF = NettoFlux(t,h,SoilPar,ModelDim,BoundaryPar)
q = Richards1DTHe(t,h,SoilPar,ModelDim,BoundaryPar);
nin = ModelDim.nin;
dzin = ModelDim.dzin;

%Calculate differential water capacity
%Cp = DiffWaterCapacity(t,h,ModelPar,ModelDim);
%SwS=1e-9;
%Cp = C+SwS;

NF = -(q(2:nin,1)-q(1:nin-1,1))./(dzin(1:nin-1,1));