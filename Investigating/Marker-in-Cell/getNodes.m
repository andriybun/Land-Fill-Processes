function modelDim = getNodes(zBottom, dz)
%
% Function generates structure containing information about geometrical
% discretization of investigated system (nodes / internodes)
%

    modelDim = struct();
    
    % Positions of internodes
    modelDim.zin = [0:-dz:zBottom]';
    % Number of internodes
    modelDim.nin = length(modelDim.zin);
    % Positions of nodes
    modelDim.zn(1, 1) = modelDim.zin(1, 1);
    modelDim.zn(2:nin-2, 1) = (modelDim.zin(2:modelDim.nin-2) + modelDim.zin(3:nin-1)) ./ 2;
    modelDim.zn(nin-1, 1) = modelDim.zin(modelDim.nin);
    % Number of nodes
    modelDim.nn = length(modelDim.zn);
    % Intervals between nodes
    modelDim.dzn = modelDim.zn(2:end, 1) - modelDim.zn(1:end-1, 1);
    % Intervals between internodes
    modelDim.dzin = modelDim.zin(2:end, 1) - modelDim.zin(1:end-1, 1);
end