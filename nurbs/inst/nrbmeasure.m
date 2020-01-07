% NRBMEASURE: Compute the distance between two given points along a NURBS curve.
% 
% Calling Sequence:
% 
%  [dist, ddistds, ddistde] = nrbmeasure (nrb)
%  [dist, ddistds, ddistde] = nrbmeasure (nrb, s, e)
%  [dist, ddistds, ddistde] = nrbmeasure (nrb, s, e, tol)
% 
% INPUT:
% 
%   nrb	: a NURBS curve, see nrbmak.
%   s   : starting point in the parametric domain.
%   e   : ending point in the parametric domain.
%   tol : tolerance for numerical quadrature, to be used in quad.
% 
% OUTPUT:
% 
%   dist   : distance between the two points along the NURBS curve.
%   ddistds: derivative of the distance function with respect to the point s.
%   ddistde: derivative of the distance function with respect to the point e.
% 
% Description:
% 
%   Compute the distance between two given points along a NURBS curve, using 
%   quad for numerical integration. The points are given by their coordinates
%   in the parametric domain.
%
% Examples:
%
%  Compute the length of a circular arc constructed as a NURBS.
%
%   c = nrbcirc (1, [0 0], 0, pi/2);
%   s = 0; e = 1;
%   l = nrbmeasure (c, s, e, 1e-7);
%
% Copyright (C) 2013 Carlo de Falco
% 
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
   
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with Octave; see the file COPYING.  If not, see
%   <http://www.gnu.org/licenses/>.  

function [dist, ddistds, ddistde] = nrbmeasure (nrb, s, e, tol)
  
  if (nargin < 4)
    tol = 1e-6;
    if (nargin < 3)
      e = 1;
      if (nargin < 2)
        s = 0;
      end
    end
  end

  nrb.knots = (nrb.knots - nrb.knots(1)) / (nrb.knots(end) - nrb.knots(1));

  if (numel (s) > 1 && isscalar (e))
    e = e * ones (size(s));
  elseif (numel (e) > 1 && isscalar (s))
    s = s * ones (size(e));
  end

  ders = nrbderiv (nrb);
  dist = arrayfun (@(x, y) quad (@(u) len (u, nrb, ders), x, ...
                                 y, tol), s, e);
  
  if (nargout > 1)
    ddistds = -len (s, nrb, ders);
    if (nargout > 2)
      ddistde = +len (e, nrb, ders);
    end
  end

end

function l = len (u, nrb, ders)
  [~, d] = nrbdeval (nrb, ders, u);
  f = d(1, :);
  g = d(2, :);
  h = d(3, :);
  l = sqrt (f.^2 + g.^2 + h.^2);
end

%!test 
%! c = nrbcirc (1, [0 0], 0,  pi/3);
%! l = nrbmeasure(c, 0, 1, 1e-7);
%! assert (l, pi/3, 1e-7)

%!test 
%! c = nrbcirc (1, [0 0], 0,  pi/2);
%! s = zeros (1, 100); e = linspace (0, 1, 100);
%! for ii = 1:100
%!   l(ii) = nrbmeasure (c, s(ii), e(ii), 1e-7);
%! endfor
%! xx = nrbeval (c, e);
%! theta = atan2 (xx(2,:), xx(1,:));
%! assert (l, theta, 1e-7)

%!test 
%! c = nrbcirc (1, [0 0], 0,  pi/2);
%! s = 0; e = linspace (0, 1, 100);
%! for ii = 1:100
%!   l(ii) = nrbmeasure (c, s, e(ii), 1e-7);
%! endfor
%! l2 = nrbmeasure (c, s, e, 1e-7);
%! assert (l, l2, eps)
