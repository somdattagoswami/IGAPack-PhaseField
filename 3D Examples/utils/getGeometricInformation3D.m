function geomInfo = getGeometricInformation3D(PHTelem, cpts, newBasisVert, p, q,r, xmin, xmax, ymin, ymax,zmin, zmax)
%calculates the "geometric information" corresponding to the new basis
%vertices

geomInfo = cell(length(newBasisVert), 1);


for indexVert=1:length(newBasisVert)
    %determine the evaluation point in the reference space
    %[-1,1]x[-1,1]x[-1,1]
    
    switch newBasisVert(indexVert)
        case 1
            uhat = 0;
            vhat = -1;
            what = 0;
        case 2
            uhat = 1;
            vhat = 0;
            what = 0;
        case 3
            uhat = 0;
            vhat = 1;
            what = 0;
        case 4
            uhat = -1;
            vhat = 0;
            what = 0;
        case 5
            uhat = 0;
            vhat = 0;
            what = -1;
        case 6
            uhat = 0;
            vhat = 0;
            what = 1;
        case 7
            uhat = 0;
            vhat = 0;
            what = 0;
        case 8
            uhat = -1;
            vhat = 0;
            what = 1;
        case 9
            uhat = -1;
            vhat = 0;
            what = -1;
        case 10
            uhat = 1;
            vhat = 0;
            what = 1;
        case 11
            uhat = 1;
            vhat = 0;
            what = -1;
        case 12
            uhat = 0;
            vhat = -1;
            what = 1;
        case 13
            uhat = 0;
            vhat = -1;
            what = -1;
        case 14
            uhat = 0;
            vhat = 1;
            what = 1;
        case 15
            uhat = 0;
            vhat = 1;
            what = -1;
        case 16
            uhat = -1;
            vhat = -1;
            what = 0;
        case 17
            uhat = 1;
            vhat = -1;
            what = 0;
        case 18
            uhat = -1;
            vhat = 1;
            what = 0;
        case 19
            uhat = 1;
            vhat = 1;
            what = 0;
            
            
    end
    
    %evaluate the Bernstein polynomials
    [B_u, dB_u] = bernstein_basis(uhat,p);
    [B_v, dB_v] = bernstein_basis(vhat,q);
    [B_w, dB_w] = bernstein_basis(what,r);
    
    Buvw = zeros((p+1)*(q+1)*(r+1),1);
    dBdu = zeros((p+1)*(q+1)*(r+1),1);
    dBdv = zeros((p+1)*(q+1)*(r+1),1);
    dBdw= zeros((p+1)*(q+1)*(r+1),1);
    ddBdudv= zeros((p+1)*(q+1)*(r+1),1);
    ddBdudw= zeros((p+1)*(q+1)*(r+1),1);
    ddBdvdw= zeros((p+1)*(q+1)*(r+1),1);
    ddBdudvdw= zeros((p+1)*(q+1)*(r+1),1);
    
    basisCounter = 0;
    for k=1:r+1
        for j=1:q+1
            for i=1:p+1
                
                basisCounter = basisCounter + 1;
                Buvw(basisCounter) = B_u(i)*B_v(j)*B_w(k);
                dBdu(basisCounter) = dB_u(i)*B_v(j)*B_w(k);
                dBdv(basisCounter) = B_u(i)*dB_v(j)*B_w(k);
                dBdw(basisCounter) = B_u(i)*B_v(j)*dB_w(k);
                ddBdudv(basisCounter) = dB_u(i)*dB_v(j)*B_w(k);
                ddBdudw(basisCounter) = dB_u(i)*B_v(j)*dB_w(k);
                ddBdvdw(basisCounter) = B_u(i)*dB_v(j)*dB_w(k);
                ddBdudvdw(basisCounter) = dB_u(i)*dB_v(j)*dB_w(k);
            end
        end
    end
    cR = PHTelem.C * Buvw;
    cdRdx = (PHTelem.C * dBdu)*2/(xmax-xmin);
    cdRdy = (PHTelem.C * dBdv)*2/(ymax-ymin);
    cdRdz = (PHTelem.C * dBdw)*2/(zmax-zmin);
    cdRdxdy = (PHTelem.C * ddBdudv)*4/((xmax-xmin)*(ymax-ymin));
    cdRdxdz = (PHTelem.C * ddBdudw)*4/((xmax-xmin)*(zmax-zmin));
    cdRdydz = (PHTelem.C * ddBdvdw)*4/((ymax-ymin)*(zmax-zmin));
    cdRdxdydz = (PHTelem.C * ddBdudvdw)*8/((xmax-xmin)*(ymax-ymin)*(zmax-zmin));
    
    coord = cR'*cpts;
    coord_u = cdRdx'*cpts;
    coord_v = cdRdy'*cpts;
    coord_w = cdRdz'*cpts;
    coord_uv = cdRdxdy'*cpts;
    coord_uw = cdRdxdz'*cpts;
    coord_vw = cdRdydz'*cpts;
    coord_uvw = cdRdxdydz'*cpts;
    
    geomInfo{indexVert} = [coord; coord_u; coord_v; coord_w; coord_uv; coord_uw; coord_vw; coord_uvw];
    
end
