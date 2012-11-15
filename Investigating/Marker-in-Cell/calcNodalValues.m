function [cNodes, sIdx] = CalcNodalValues(zMark, cMark, ModelDim)
% Function:
%   Calc nodal values
%
% Purpose:
%   Function to calculate nodal values given values at a set of markers.
%   It uses markers with coordinates within half grid step from a given
%   node (between nearest internodes) to calculate values at any node.
%
% Parameters:
%   zMark - vector of coordinates of markers.
%   cMark - vector of concentrations at markers. Must be of the same
%           size as zMark.
%   ModelDim - structure containing all information about the
%              geometrical properties of a system. It must contain following
%   members:
%       - zn (coordinates of nodes);
%       - zin (coordinates of internodes);
%       - nn (number of nodes);
%
% Return:
%   cNodes - concentration of values at given nodes

    % Sort markers
    [zMarkSorted, sIdx] = sort(-zMark);
    cMarkSorted = cMark(sIdx);

    % Boundaries of cells to consider
    zInDown = -ModelDim.zin(1:end-1);
    zInUp = -ModelDim.zin(2:end);
    
    % Nodal concentrations
    cNodes = zeros(ModelDim.nn, 1);
    
    for idx = 1:ModelDim.nn
        % Markers that are located near current node
        conIdx = (zMarkSorted >= zInDown(idx)) & (zMarkSorted <= zInUp(idx));
        if any(conIdx)
            % Distance from marker to current node
            dzMark = abs(zMarkSorted(conIdx) + ModelDim.zn(idx));
            % Maximum distance from marker to current node
            dzMarkMax = max(abs(-ModelDim.zn(idx) - zInUp(idx)), abs(-ModelDim.zn(idx) - zInDown(idx)));
            % Weights of markers for current node
            wM = 1 - dzMark ./ dzMarkMax;
            % Concentration at the node
            cNodes(idx) = sum(wM .* cMarkSorted(conIdx)) / sum(wM);
        else
            cNodes(idx) = 0;
        end
    end
    
end