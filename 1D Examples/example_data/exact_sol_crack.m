function [y,dy,f] = exact_sol_crack(x)
% Defines the exact solution from equation (63) in Schillinger - IGA Phase
% field collocation paper

xneg = (x<0);
xpos = (x>=0);

y = zeros(size(x));
y(xneg) = 1/pi^2*sin(pi*x(xneg))-(1+x(xneg))/pi;
y(xpos) = 1/pi^2*sin(pi*x(xpos))+(1-x(xpos))/pi;

dy = 1/pi*cos(pi*x)-1/pi;
ddy = -sin(pi*x);

f = -ddy;
end
