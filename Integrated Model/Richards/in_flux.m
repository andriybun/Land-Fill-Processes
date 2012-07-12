function in_flx = in_flux(t)
    in_flx = 0;
    if (t == 2)
        in_flx = 1e-3;
    end
%     in_flx = 0;
%     shift = 50;
%     if (t >= shift) && (t < shift + 50)
%         in_flx = (t - shift) / 50 * 9e-3;
%     end
%     if (t >= shift + 50) && (t < shift + 100)
%         in_flx = 9e-3;
%     end
% %     if (t >= shift + 100) && (t < shift + 800)
% %         in_flx = 2.5e-4;
% %     end
end
