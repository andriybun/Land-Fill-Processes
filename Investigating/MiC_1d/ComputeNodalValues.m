function [cNodes, sIdx] = ComputeNodalValues(zMark, cMark, ModelDim)
% Function:
%   Compute nodal values
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
%   zNode - coordinates of nodes.
%
% Return:
%   cNodes - concentration of values at given nodes.
%   sIdx   - indices of elements in zMark while sorted.

    INTERPOLATION_METHOD = 'pchip';

    % Sort markers
    [zMarkSorted, sIdx] = sort(-zMark);

    %% Alternative calculation using interp1:
%     cNodes = interp1(zMark(sIdx), cMark(sIdx), ModelDim.zn, INTERPOLATION_METHOD, 'extrap');
    
    %% Alternative calculation using gridfit:
%     xAux = 1:numel(zMark);
%     cNodes = gridfit(-zMark(sIdx), xAux, cMark(sIdx), -ModelDim.zn, xAux);

    %% Alternative using average
    cMarkSorted = cMark(sIdx);

    % Boundaries of cells to consider
    zInDown = -ModelDim.zin(1:end-1);
    zInUp = -ModelDim.zin(2:end);
    
    % Nodal concentrations
    nNodes = numel(ModelDim.zin) - 1;
    cNodes = zeros(nNodes, 1);
    
    for idx = 1:nNodes
        % Markers that are located near current node
        isInNode = (zMarkSorted >= zInDown(idx)) & (zMarkSorted <= zInUp(idx));
        if any(isInNode)
            % Concentration at the node
            cNodes(idx) = mean(cMarkSorted(isInNode));
        else
            cNodes(idx) = 0;
        end
    end

end
