function vMark = ComputeMarkerValues(zNode, vNode, zMark)
% Function:
%   Compute marker values
%
% Purpose:
%   Function to calculate nodal values given values at a set of markers.
%   It uses markers with coordinates within half grid step from a given
%   node (between nearest internodes) to calculate values at any node.
%
% Parameters:
%   zNode - coordinates of nodes.
%   vNodes - values at given nodes. Must be of the same size as modelDim.zn.
%   zMark - vector of coordinates of markers.
%
% Return:
%   vMark - vector of values at markers' coordinates.

    INTERPOLATION_METHOD = 'linear';
    vMark = interp1(zNode, vNode, zMark, INTERPOLATION_METHOD, 'extrap');
    
end