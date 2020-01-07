function S = basiskntins (deg,t,u)

% Compute the coefficient matrix for non-uniform B-splines subdivision.
%
% This represents the B-spline basis given by a coarse knot vector
%  in terms of the B-spline basis of a finer knot vector.
%
% The function is implemented for the univariate case. It is based on 
%  the paper:
%
% G. Casciola, L. Romani, A general matrix representation for non-uniform
%  B-spline subdivision with boundary control, ALMA-DL, University of Bologna (2007)
%
% Calling Sequence:
% 
%    S = basiskntins (deg, t, u);
%
%    INPUT:
%   
%      deg - degree of the first knot vector
%      t   - coarse knot vector
%      u   - fine knot vector
%   
%    OUTPUT:
%   
%      B - Value of the basis functions at the points
%          size(B)=[numel(u),(p+1)] for curves
%          or [numel(u)*numel(v), (p+1)*(q+1)] for surfaces
%
%      N - Indices of the basis functions that are nonvanishing at each
%          point. size(N) == size(B)
%   
%    Copyright (C) 2015 Rafael Vazquez
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

nt = length(t);
nu = length(u);
S = sparse (nu-deg-1,nt-deg-1);
[t_mult,t_single,nt_s] = knot_mult(deg,t);
[u_mult,u_single,nu_s] = knot_mult(deg,u);
st = deg+1;
su = deg+1;
row = 1;
col = 1;
Sl = bs2bs(deg,t,u,st,su);
S(row:deg+row,col:deg+col) = Sl;
t_single(nt+1) = t(nt-deg);
i = 1;

for j=1:nu_s
  if (u_single(j) == t_single(i))
    st = st+t_mult(i);
    col = col+t_mult(i);
    i = i+1;
  end
  su = su+u_mult(j);
  row = row+u_mult(j);
  Sl = bs2bs(deg,t,u,st,su);
  S(row:deg+row,col:deg+col) = Sl;
end

end


function [t_mult,t_single,nt_s] = knot_mult(d,t)
epsilon = 1e-14 * (t(end) - t(1));

nt = length(t);
nt_s = 0;
m = 1;
for i = d+2:nt-d-1
  if ((t(i+1) - t(i)) > epsilon)
    nt_s = nt_s+1;
    t_mult(nt_s) = m;
    t_single(nt_s) = t(i);
    m=1;
  else
    m = m+1;
  end
end
t_single(nt_s+1)=t(nt-d);
t_mult(nt_s+1)=0;

end

function S = bs2bs(d,t,u,k,l)

S = zeros(d+1);
S(1,:) = bs2bs_first_row(d,t,u,k,l);

for ir=1:d
  S(ir+1,:) = bs2bs_i_row(d,t,u,k,l,ir,S(ir,:));
end
end

function S = bs2bs_first_row(d,t,u,k,l)

S = eye(1,d+1);
for h=1:d
  beta_2=0.0;
  uu=u(l+1-h);
  for j=h:-1:1
    d1=uu-t(k+j-h);
    d2=t(k+j)-uu;
    beta_1=S(j)/(d2+d1);
    S(j+1)=d1*beta_1+beta_2;
    beta_2=d2*beta_1;
  end
  S(1)=beta_2;
end
end

function Si = bs2bs_i_row(d,t,u,k,l,ir,S)

Si(1) = S(1)*(t(k+1)-u(l+ir))/(t(k+1)-u(l+ir-d));

for j=1:d
  den=t(k+j+1)-u(l+ir-d);
  fact=(t(k+j+1)-t(k+j-d))/(t(k+j)-t(k+j-d-1));
  Si(j+1)=(fact*(S(j)*(u(l+ir)-t(k+j-d-1))-Si(j) * ...
    (u(l+ir-d)-t(k+j-d-1)))+S(j+1)*(t(k+j+1)-u(l+ir)))/den;
end
end

%!test
%! knt1 = [0 0 0 1/2 1 1 1];
%! knt2 = [0 0 0 1/4 1/2 3/4 1 1 1];
%! C = basiskntins (2, knt1, knt2);
%! assert (full(C), [1 0 0 0; 1/2 1/2 0 0; 0 3/4 1/4 0; 0 1/4 3/4 0; 0 0 1/2 1/2; 0 0 0 1]);

%!test
%! crv = nrbtestcrv;
%! crv2 = nrbkntins (crv, [0.1, 0.3, 0.4, 0.5, 0.6, 0.8, 0.96 0.98]);
%! C = basiskntins (crv.order-1,crv.knots,crv2.knots);
%! for ii = 1:4
%!   assert (max (abs(C*crv.coefs(ii,:)' - crv2.coefs(ii,:)')) < 1e-14 )
%! end

%!test
%! crv = nrbtestcrv;
%! crv2 = nrbkntins (crv, [0.50000001, 0.5000001, 0.500001, 0.50001, 0.5001]);
%! C = basiskntins (crv.order-1,crv.knots,crv2.knots);
%! for ii = 1:4
%!   assert (max (abs(C*crv.coefs(ii,:)' - crv2.coefs(ii,:)')) < 1e-14 )
%! end

%!test
%! crv = nrbtestcrv;
%! crv2 = nrbkntins (crv, [0.1, 0.3, 0.4, 0.5, 0.6, 0.8, 0.96 0.98]);
%! C = basiskntins (crv.order-1,crv.knots,crv2.knots);
%! x = linspace (0, 1, 10);
%! s = findspan (crv.number-1, crv.order-1, x, crv.knots);
%! s2 = findspan (crv2.number-1, crv2.order-1, x, crv2.knots);
%! N = basisfun (s, x, crv.order-1, crv.knots);
%! N2 = basisfun (s2, x, crv2.order-1, crv2.knots);
%! c = numbasisfun (s, x, crv.order-1, crv.knots) + 1;
%! c2 = numbasisfun (s2, x, crv2.order-1, crv2.knots) + 1;
%! for ii = 1:numel(x)
%!   assert (abs(N2(ii,:) * C(c2(ii,:),c(ii,:)) - N(ii,:)) < 1e-14)
%! end
