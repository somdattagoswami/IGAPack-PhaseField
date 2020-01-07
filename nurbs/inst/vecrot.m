function rx = vecrot(angle, vector)
% 
% VECROT: Transformation matrix for a rotation around the axis given by a vector. 
% 
% Calling Sequence:
% 
%   rx = vecrot (angle, vector);
% 
% INPUT:
% 
%   angle : rotation angle defined in radians
%   vector: vector defining the rotation axis
%
% OUTPUT:
% 
%   rx: (4x4) Transformation matrix.
% 
% 
% Description:
% 
%   Return the (4x4) Transformation matrix for a rotation about the axis
%   defined by vector, and by the given angle.
% 
% See also:
% 
%    nrbtform
%
%    Copyright (C) 2011 Rafael Vazquez
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

% Normalize the vector
vec = vector / norm (vector);

sn = sin (angle);
cn = cos (angle);
rx = [cn+vec(1)^2*(1-cn), vec(1)*vec(2)*(1-cn)-vec(3)*sn, vec(1)*vec(3)*(1-cn)+vec(2)*sn, 0; 
    vec(1)*vec(2)*(1-cn)+vec(3)*sn, cn+vec(2)^2*(1-cn), vec(2)*vec(3)*(1-cn)-vec(1)*sn, 0; 
    vec(1)*vec(3)*(1-cn)-vec(2)*sn, vec(2)*vec(3)*(1-cn)+vec(1)*sn, cn+vec(3)^2*(1-cn), 0; 
    0 0 0 1];

end
