function q = Richards1DTHe(t,h,SoilPar,ModelDim,BoundaryPar)
%Function for calculating rates for a Heat Flow in a spatial domain
%qH are the heat fluxes at all internodal points in the spatial domain;
%t is time
%T are the temperatures in all cells in the spatial domain
%ModelPar contains the Model Parameters (thermal diffusivity of all
%InterNodal points)
%ModelDim contains the ModelDiscretization (coordinates of nodes and 
% internodes, deltas of nodal points and
%internodal points)

nn = ModelDim.nn;
nin = ModelDim.nin;
dzn = ModelDim.dzn;
dzin = ModelDim.dzin;

%Calculate Derived States (total head etc.)
hin(1,1) = h(1,1);
hin(2:nn,1) = (dzin(1:nn-1,1).*h(1:nn-1,1) + dzin(2:nn,1).*h(2:nn,1))./...
    (dzin(1:nn-1,1)+dzin(2:nn));
hin(nn+1,1) = h(nn,1);
K = HyC_function(hin,SoilPar,ModelDim);

%upper boundary zero flux
q(1,1) = BoundaryPar.qTop(t);

q(2:nin-1,1)=-K(2:nin-1,1).*((h(2:nn,1)-h(1:nn-1,1))./dzn(1:nn-1,1) + 1);

% Lower/Right boundary conditions
% unit gradient (gravity drainage (zero pressure gradient) ...

q(nin,1)= BoundaryPar.Ksurf.*(BoundaryPar.hamb-h(nn,1));

