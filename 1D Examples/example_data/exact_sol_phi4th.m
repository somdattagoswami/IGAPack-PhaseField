function [y] = exact_sol_phi4th(x,Fract)

y = zeros(size(x));
y = exp(-abs(x-0.0)/Fract.constl).*(1+abs(x-0.0)/Fract.constl);

end