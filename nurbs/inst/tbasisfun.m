function [N, Nder] = tbasisfun (u, p, U)
%
% TBASISFUN: Compute a B- or T-Spline basis function, and its derivatives, from its local knot vector.
%
% usage:
%
% [N, Nder] = tbasisfun (u, p, U)
% [N, Nder] = tbasisfun ([u; v], [p q], {U, V})
% [N, Nder] = tbasisfun ([u; v; w], [p q r], {U, V, W})
% 
% INPUT:
%
%  u or [u; v] : points in parameter space where the basis function is to be
%  evaluated 
%  
%  U or {U, V} : local knot vector
%
% p or [p q] : polynomial order of the basis function
%
% OUTPUT:
%
%  N    : basis function evaluated at the given parametric points 
%  Nder : basis function gradient evaluated at the given parametric points 
%
%    Copyright (C) 2009 Carlo de Falco
%    Copyright (C) 2012 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
  if (~ iscell (U))
    U = sort (U);
    if (numel (U) ~= p+2)
      error ('tbasisfun: knot vector and degree do not correspond')
    end
    
    if (nargout == 1)
      N = onebasisfun__ (u, p, U);
    else
      [N, Nder] = onebasisfunder__ (u, p, U);
    end
    
  elseif (size(U,2) == 2)
    U{1} = sort(U{1}); U{2} = sort(U{2});
    if (numel(U{1}) ~= p(1)+2 || numel(U{2}) ~= p(2)+2)
      error ('tbasisfun: knot vector and degree do not correspond')
    end
    
    if (nargout == 1)
      Nu = onebasisfun__ (u(1,:), p(1), U{1});
      Nv = onebasisfun__ (u(2,:), p(2), U{2});

      N = Nu.*Nv;
    elseif (nargout == 2)
      [Nu, Ndu] = onebasisfunder__ (u(1,:), p(1), U{1});
      [Nv, Ndv] = onebasisfunder__ (u(2,:), p(2), U{2});

      N = Nu.*Nv;
      Nder(1,:) = Ndu.*Nv;
      Nder(2,:) = Nu.*Ndv;
    end
    
  elseif (size(U,2) == 3)
    U{1} = sort(U{1}); U{2} = sort(U{2}); U{3} = sort(U{3});
    if (numel(U{1}) ~= p(1)+2 || numel(U{2}) ~= p(2)+2 || numel(U{3}) ~= p(3)+2)
      error ('tbasisfun: knot vector and degree do not correspond')
    end

    if (nargout == 1)
      Nu = onebasisfun__ (u(1,:), p(1), U{1});
      Nv = onebasisfun__ (u(2,:), p(2), U{2});
      Nw = onebasisfun__ (u(3,:), p(3), U{3});

      N = Nu.*Nv.*Nw;
    else
      [Nu, Ndu] = onebasisfunder__ (u(1,:), p(1), U{1});
      [Nv, Ndv] = onebasisfunder__ (u(2,:), p(2), U{2});
      [Nw, Ndw] = onebasisfunder__ (u(3,:), p(3), U{3});
      
      N = Nu.*Nv.*Nw;
      Nder(1,:) = Ndu.*Nv.*Nw;
      Nder(2,:) = Nu.*Ndv.*Nw;
      Nder(3,:) = Nu.*Nv.*Ndw;
    end
  end

end

%!demo
%! U = {[0 0 1/2 1 1], [0 0 0 1 1]};
%! p = [3, 3];
%! [X, Y] = meshgrid (linspace(0, 1, 30));
%! u = [X(:), Y(:)]';
%! N = tbasisfun (u, p, U);
%! surf (X, Y, reshape (N, size(X)))
%! title('Basis function associated to a local knot vector')
%! hold off

%!test
%! U = [0 1/2 1];
%! p = 1;
%! u = [0.3 0.4 0.6 0.7];
%! [N, Nder] = tbasisfun (u, p, U);
%! assert (N, [0.6 0.8 0.8 0.6], 1e-12);
%! assert (Nder, [2 2 -2 -2], 1e-12);

%!test
%! U = {[0 1/2 1] [0 1/2 1]};
%! p = [1 1];
%! u = [0.3 0.4 0.6 0.7; 0.3 0.4 0.6 0.7];
%! [N, Nder] = tbasisfun (u, p, U);
%! assert (N, [0.36 0.64 0.64 0.36], 1e-12);
%! assert (Nder, [1.2 1.6 -1.6 -1.2; 1.2 1.6 -1.6 -1.2], 1e-12);

%!test
%! U = {[0 1/2 1] [0 1/2 1] [0 1/2 1]};
%! p = [1 1 1];
%! u = [0.4 0.4 0.6 0.6; 0.4 0.4 0.6 0.6; 0.4 0.6 0.4 0.6];
%! [N, Nder] = tbasisfun (u, p, U);
%! assert (N, [0.512 0.512 0.512 0.512], 1e-12);
%! assert (Nder, [1.28 1.28 -1.28 -1.28; 1.28 1.28 -1.28 -1.28; 1.28 -1.28 1.28 -1.28], 1e-12);
