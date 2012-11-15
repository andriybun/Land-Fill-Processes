close all
x = 0:0.1:10;

pH1 = 9;
pH2 = 6;

xdiff = (pH1-pH2)/2;
xscale = 2/xdiff;
shift = xdiff+pH2;

y =  1/2 * erfc(-xscale * (x - shift));

plot(x,y)

