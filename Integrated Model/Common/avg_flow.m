% This is a function to calculate average flow for a sequence of time
% intervals with a given flow at the intervals' boundaries.
function res = avg_flow(flow)
    d = ndims(flow);
    if d == 2
        % 1D flow
        flow = cat(1, flow, flow(end, :));
        res = (flow(2:end, :) + flow(1:end-1, :)) ./ 2;
    elseif d == 3
        % 2D or 3D flow
        flow = cat(1, flow, flow(end, :, :));
        res = (flow(2:end, :, :) + flow(1:end-1, :, :)) ./ 2;
    end
%     plot (flow);
%     hold on;
%     plot (res, 'r');
%     hold off;
end