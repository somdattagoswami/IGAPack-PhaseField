function geomInfo = getGeometricInformation1D(PHTelem,cpts,p,umin,umax)
%calculates the "geometric information" corresponding to the new basis
%vertices

pcpts = zeros(size(cpts));
pcpts(:,1) = cpts(:,1).*cpts(:,2);
pcpts(:,2) = cpts(:,2);
pcpts(:,3:end) = cpts(:,3:end);

uhat = 0;
%evaluate the Bernstein polynomials
[B_u,dB_u] = bernstein_basis(uhat,p);
cR = PHTelem.C * B_u';
cdRdx = (PHTelem.C * dB_u')*2/(umax-umin);

coord = cR'*pcpts;
coord_u = cdRdx'*pcpts;
geomInfo = [coord; coord_u];

end
