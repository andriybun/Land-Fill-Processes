function y = smoothstep(pH1,pH2,x)

xx = zeros(length(pH1),1);
xx(:) = x;
xdiff = (pH1-pH2)./2;
xscale = 2./xdiff;
shift = xdiff+pH2;

y =  1/2 .* erfc(-xscale .* (xx - shift));
end

