function tvol = nrbpermute (vol, ord)
% 
% NRBPERMUTE: Rearrange the directions of a NURBS volume or surface.
% 
% Calling Sequence:
% 
%   tvol = nrbpermute(vol,order)
%
% INPUT:
% 
%   vol		: NURBS volume or surface, see nrbmak.
%   order   : the order to rearrange the directions of the NURBS entity.
%
% OUTPUT:
% 
%   tvol	: NURBS volume or surface with rearranged directions.
% 
% Description:
% 
%   Utility function that rearranges the directions of a NURBS volume or
%   surface. For surfaces, nrbpermute(srf,[1 2]) is the same as
%   nrbtransp(srf). NURBS curves cannot be rearranged.
%
% Example:
%
%    nrbpermute (vol, [1 3 2])
%
%    Copyright (C) 2013 Rafael Vazquez
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

if (~iscell(vol.knots))
  error('A NURBS curve cannot be rearranged.');
end

tvol = nrbmak (permute (vol.coefs, [1, ord+1]), {vol.knots{ord}});

%!demo
%! vol = nrbrevolve (nrb4surf ([1 0], [2 0], [1 1], [2 1]), [0 0 0], [0 1 0], pi/8);
%! nrbplot(vol,[5 10 20]);
%! title('NURBS volume and the same after reordering the directions')
%! hold on
%! vol.coefs(1,:,:) = vol.coefs(1,:,:) + 2;
%! vol = nrbpermute(vol,[2 3 1]);
%! nrbplot(vol,[5 10 20]);
%! hold off

%!test
%! vol = nrbrevolve (nrb4surf ([1 0], [2 0], [1 1], [2 1]), [0 0 0], [0 1 0], pi/8);
%! perm1 = [1 3 2];
%! perm2 = [2 1 3];
%! vol2 = nrbpermute (vol, perm1);
%! vol3 = nrbpermute (vol, perm2);
%! assert (vol.number(perm1), vol2.number)
%! assert (vol.order(perm1), vol2.order)
%! assert ({vol.knots{perm1}}, vol2.knots)
%! assert (permute(vol.coefs, [1, perm1+1]), vol2.coefs)
%! assert (vol.number(perm2), vol3.number)
%! assert (vol.order(perm2), vol3.order)
%! assert ({vol.knots{perm2}}, vol3.knots)
%! assert (permute(vol.coefs, [1, perm2+1]), vol3.coefs)
