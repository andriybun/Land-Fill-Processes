function ModelDim = InitializeNodes(varName, varBottom, dVar)
%
% Function generates structure containing information about geometrical
% discretization of investigated system (nodes / internodes)
% Function generates a list of field names corresponding to a variable
% characterizing one spatial dimension.
% Fields generated are:
%   - positions of internodes (from 0 to varBottom with a step dVar)
%   - positions of nodes
%   - intervals between internodes
%   - intervals between nodes
%   - number of internodes
%   - number of nodes
%

    % Templates for field names. Field names are generated using passed
    % varName.
    fnTemplates = {'%sin', '%sn', 'd%sin', 'd%sn', '%snin', '%snn',};
    vin  = sprintf(fnTemplates{1}, varName);
    vn   = sprintf(fnTemplates{2}, varName);
    dvin = sprintf(fnTemplates{3}, varName);
    dvn  = sprintf(fnTemplates{4}, varName);
    vnin = sprintf(fnTemplates{5}, varName);
    vnn  = sprintf(fnTemplates{6}, varName);
    
    % Initialize output structure
    ModelDim = struct();
    
    % Initialize its fields
    % Positions of internodes (equally spaced on interval [0, varBottom])
    ModelDim.(vin) = [0:-dVar:varBottom]';
    % Number of internodes
    ModelDim.(vnin) = length(ModelDim.(vin));
    % Positions of nodes (in between internodes, first and last nodes are
    % aligned with first and last internodes respectively)
    ModelDim.(vn)(1, 1) = ModelDim.(vin)(1, 1);
    ModelDim.(vn)(2:ModelDim.(vnin)-2, 1) = (ModelDim.(vin)(2:ModelDim.(vnin)-2) + ModelDim.(vin)(3:ModelDim.(vnin)-1)) ./ 2;
    ModelDim.(vn)(ModelDim.(vnin)-1, 1) = ModelDim.(vin)(ModelDim.(vnin));
    % Number of nodes
    ModelDim.(vnn) = length(ModelDim.(vn));
    % Intervals between nodes
    ModelDim.(dvn) = ModelDim.(vn)(2:end, 1) - ModelDim.(vn)(1:end-1, 1);
    % Intervals between internodes
    ModelDim.(dvin) = ModelDim.(vin)(2:end, 1) - ModelDim.(vin)(1:end-1, 1);
end