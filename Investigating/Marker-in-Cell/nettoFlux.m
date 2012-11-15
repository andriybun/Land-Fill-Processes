function nf = nettoFlux(t, C, SoilPar, ModelDim, BoundaryPar)
    qH = soluteFlow1d(t, C, SoilPar, ModelDim, BoundaryPar);
    nin = ModelDim.nin;
    dzin = ModelDim.dzin;

    nf = -(qH(2:nin, 1) - qH(1:nin-1, 1)) ./ (dzin(1:nin-1));
end