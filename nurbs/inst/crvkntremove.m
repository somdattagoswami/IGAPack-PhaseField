

function [rcrv, t] = crvkntremove (crv, u, r, s, num, d)
% 
% CRVKNTREMOVE: Remove one knot from the knot-vector of a NURBS curve.
% 
% Calling Sequence:
% 
%   [rcrv, remflag] = crvkntremove (crv, u, r, s, num, d);
% 
% INPUT:
% 
%   crv		: NURBS curve, see nrbmak.
% 
%   u           : knot to be removed.
% 
%   r           : index of the knot to be removed.
% 
%   s           : multiplicity of the knot to be removed.
% 
%   num         : number of knot removals requested.
%
%   d           : curve deviation tolerance.
%
% OUTPUT:
%
%   rcrv        : new NURBS structure for the curve with knot u remuved.
% 
%   t           : actual number of knot removals performed.
% 
% 
% 
% DESCRIPTION:
% 
%   Remove knot u from the NURBS curve crv at most num times. 
%   Check that the maximum deviation of the curve be less than d.
%   Based on algorithm A5.8 NURBS Book (pag183)
%  
% SEE ALSO:
% 
%   nrbkntins
%
%    Copyright (C) 2013 Jacopo Corno
%    Copyright (C) 2013 Carlo de Falco
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

  [U, Pw, t] = RemoveCurveKnot (crv.number, crv.order - 1, crv.knots, ...
                                crv.coefs, u, r, s, num, d);
  rcrv = nrbmak (Pw, U);

end

%!test 
%! crv  = nrbdegelev (nrbline (), 3);
%! acrv = nrbkntins (crv, [.11 .11 .11]);
%! [rcrv, t] = crvkntremove (acrv, .11, 8, 3, 3, 1e-10);
%! assert (crv.knots, rcrv.knots, 1e-10);
%! assert (t, 3);

%!test 
%! crv  = nrbcirc ();
%! acrv = nrbkntins (crv, [.3 .3]);
%! [rcrv, t] = crvkntremove (acrv, .3, 7, 2, 2, 1e-10);
%! assert (crv.knots, rcrv.knots, 1e-10);
%! assert (t, 2);

function [U, Pw, t] = RemoveCurveKnot (n, p, U, Pw, u, r, s, num, d)
% see algorithm A5.8 NURBS Book (pag183)
  
  w = min (Pw(4,:));
  Pmax = max (sqrt (sum (Pw.^2, 1)));
  TOL = d*w / (1 + Pmax);
  
  m     = n + p + 1;
  ord   = p + 1;
  fout  = (2*r - s - p) / 2; % first control point out
  last  = r - s;
  first = r - p;
  
  temp = zeros (4, 2*p + 1);
  
  for t = 0:num-1
    off = first - 1; % diff in index between temp and P
    temp(:,1) = Pw(:,off);
    temp(:,last+1-off+1) = Pw(:,last+1);
    i   = first;
    j   = last;
    ii  = 1;
    jj  = last - off;
    remflag = 0;
    while (j - i > t)
      % compute new control points for one removal step
      alfi = (u-U(i)) / (U(i+ord+t)-U(i));
      alfj = (u-U(j-t)) / (U(j+ord)-U(j-t));
      temp(:,ii+1) = (Pw(:,i)-(1.0-alfi).*temp(:,ii-1+1))./alfi;
      temp(:,jj+1) = (Pw(:,j)-alfj.*temp(:,jj+1+1))./(1.0-alfj);
      i   = i + 1;
      ii  = ii + 1;
      j   = j - 1;
      jj  = jj - 1;
    end
    if (j - i <= t)
      % check if knot removable
      if (norm (temp(:,ii-1+1) - temp(:,jj+1+1)) <= TOL)
        remflag = 1;
      else
        alfi = (u-U(i)) / (U(i+ord+t)-U(i));
        if (norm (Pw(:,i) - (alfi.*temp(:,ii+t+1+1) + ...
                           (1-alfi).*temp(:,ii-1+1))) <= TOL)
          remflag = 1;
        end%if
      end%if  
    end%if
    if (remflag == 0)
      break; % cannot remove any more knots -> get out of for loop
    else
      % successful removal -> save new control points
      i = first;
      j = last;
      while (j - i > t)
        Pw(:,i) = temp(:,i-off+1);
        Pw(:,j) = temp(:,j-off+1);
        i = i + 1;
        j = j - 1;
      end
    end%if
    first = first - 1;
    last = last + 1;
    t = t + 1;
  end % end of for loop

  if (t == 0)
    return;
  end%if
  
  % shift knots
  for k = r+1:m
    U(k-t) = U(k);
  end
  U = U(1:end-t);
  
  j = floor(fout);
  i = j;
  for k = 1:t-1
    if (mod (k, 2) == 1)
      i = i+1;
    else
      j = j-1;
    end%if
  end
  
  % shift points
  for k = i+1:n
    Pw(:,j) = Pw(:,k);
    j = j+1;
  end
  Pw = Pw(:,1:end-t);
  
  return;
end

