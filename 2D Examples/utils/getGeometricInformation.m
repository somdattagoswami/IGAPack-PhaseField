function geomInfo = getGeometricInformation(PHTelem, cpts, newBasisVert, p, q, xmin, xmax, ymin, ymax)
%calculates the "geometric information" corresponding to the new basis
%vertices

geomInfo = cell(length(newBasisVert), 1);
pcpts = zeros(size(cpts));
pcpts(:,1) = cpts(:,1).*cpts(:,3);
pcpts(:,2) = cpts(:,2).*cpts(:,3);
pcpts(:,3) = cpts(:,3);
pcpts(:,4:end) = cpts(:,4:end);

for indexVert=1:length(newBasisVert)
    %determine the evaluation point in the reference space [-1,1]x[-1,1]
   
    switch newBasisVert(indexVert)
        case 1
            uhat = 0;
            vhat = -1;
        case 2
            uhat = 1;
            vhat = 0;
        case 3
            uhat = 0;
            vhat = 1;
        case 4
            uhat = -1;
            vhat = 0;
        case 5
            uhat = 0;
            vhat = 0;
    end
    
    %evaluate the Bernstein polynomials
    [B_u, dB_u] = bernstein_basis(uhat,p);
    [B_v, dB_v] = bernstein_basis(vhat,q);
    
    Buv = zeros((p+1)*(q+1),1);
    dBdu = zeros((p+1)*(q+1),1);
    dBdv = zeros((p+1)*(q+1),1);
    ddBdudv = zeros((p+1)*(q+1),1);
    
    basisCounter = 0;
    for j=1:q+1
        for i=1:p+1
            basisCounter = basisCounter + 1;
            Buv(basisCounter) = B_u(i)*B_v(j);
            dBdu(basisCounter) = dB_u(i)*B_v(j);
            dBdv(basisCounter) = B_u(i)*dB_v(j);
            ddBdudv(basisCounter) = dB_u(i)*dB_v(j);
        end
    end   
    
    cR = PHTelem.C * Buv;  
    cdRdx = (PHTelem.C * dBdu)*2/(xmax-xmin);
    cdRdy = (PHTelem.C * dBdv)*2/(ymax-ymin);
    cdRdxdy = (PHTelem.C * ddBdudv)*4/((xmax-xmin)*(ymax-ymin));
      
    coord = cR'*pcpts;
    coord_u = cdRdx'*pcpts;
    coord_v = cdRdy'*pcpts;
    coord_uv = cdRdxdy'*pcpts;
    geomInfo{indexVert} = [coord; coord_u; coord_v; coord_uv];
  
end
