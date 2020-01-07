function plot3D(PHTelem,tdisp,sizeBasis,geometry,controlPts,vtuFile,fudge)
% Evaluate the displacements, Strain and Potential values at the element vertex

if nargin < 7
   fudge = 0.0;
end

nx = 2;
ny = 2;
nz = 2;
p = geometry.p;
q = geometry.q;
r = geometry.r;
px = linspace(-1+fudge,1-fudge,nx);
py = linspace(-1+fudge,1-fudge,ny);
pz = linspace(-1+fudge,1-fudge,nz);
noGpEle = nx*ny*nz;% Number of Gauss points per element

noElems = geometry.numberElements;
physcoord = zeros(3,noGpEle,noElems);  % Global coords of Gauss points
dispcoord     = zeros(3,noGpEle,noElems);  % Displacements of Gauss points
epsiloncoord = zeros(6,noGpEle,noElems);% Strain of Gauss points
wcoord = zeros(1,noGpEle,noElems);  % Phase Field at Gauss points

% 1D Bernstein polynomials evaluated at the Gauss points on the master element
[B_u,dB_u] = bernstein_basis(px,p);
[B_v,dB_v] = bernstein_basis(py,q);
[B_w,dB_w] = bernstein_basis(pz,r);

dBdu = zeros(nx,ny,nz,(p+1)*(q+1)*(r+1));
dBdv = zeros(nx,ny,nz,(p+1)*(q+1)*(r+1));
dBdw = zeros(nx,ny,nz,(p+1)*(q+1)*(r+1));
R = zeros(nx,ny,nz,(p+1)*(q+1)*(r+1));

basisCounter = 0;
for k=1:r+1
    for j=1:q+1
        for i=1:p+1
            basisCounter = basisCounter + 1;
            for kk=1:nz
                for jj=1:ny
                    for ii=1:nx
                        R(ii,jj,kk,basisCounter) = B_u(ii,i)*B_v(jj,j)*B_w(kk,k);
                        dBdu(ii,jj,kk,basisCounter) = dB_u(ii,i)*B_v(jj,j)*B_w(kk,k);
                        dBdv(ii,jj,kk,basisCounter) = B_u(ii,i)*dB_v(jj,j)*B_w(kk,k);
                        dBdw(ii,jj,kk,basisCounter) = B_u(ii,i)*B_v(jj,j)*dB_w(kk,k);
                    end
                end
            end
        end
    end
end

elementCounter = 1;
for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)            

            ximin = PHTelem{indexPatch}(i).vertex(1);
            ximax = PHTelem{indexPatch}(i).vertex(4);
            etamin = PHTelem{indexPatch}(i).vertex(2);
            etamax = PHTelem{indexPatch}(i).vertex(5);
            zetamin = PHTelem{indexPatch}(i).vertex(3);
            zetamax = PHTelem{indexPatch}(i).vertex(6);
            
            nument = size(PHTelem{indexPatch}(i).C,1);
            nodes = PHTelem{indexPatch}(i).nodes(1:nument);
            sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            
            cpts = controlPts{indexPatch}(nodes,1:3);
            wgts = controlPts{indexPatch}(nodes,4);
            
            gp = 1;
            for kk=1:nz
                for jj=1:ny
                    for ii=1:nx
                        dRdx = (PHTelem{indexPatch}(i).C)*squeeze(dBdu(ii,jj,kk,:))*2/(ximax-ximin);
                        dRdy = (PHTelem{indexPatch}(i).C)*squeeze(dBdv(ii,jj,kk,:))*2/(etamax-etamin);
                        dRdz = (PHTelem{indexPatch}(i).C)*squeeze(dBdw(ii,jj,kk,:))*2/(zetamax-zetamin);
                        RR = (PHTelem{indexPatch}(i).C)*squeeze(R(ii,jj,kk,:));
                        RR = RR .* wgts;
                        dRdx = dRdx .* wgts;
                        dRdy = dRdy .* wgts;
                        dRdz = dRdz .* wgts;
                        
                        w_sum = sum(RR);
                        dw_xi = sum(dRdx);
                        dw_eta = sum(dRdy);
                        dw_zeta = sum(dRdz);
                        
                        dRdx = dRdx/w_sum - RR*dw_xi/w_sum^2;
                        dRdy = dRdy/w_sum - RR*dw_eta/w_sum^2;
                        dRdz = dRdz/w_sum - RR*dw_zeta/w_sum^2;
                        RR = RR/w_sum;
                        
                        dR  = [dRdx';dRdy';dRdz'];
                        dxdxi = dR*cpts;
                        coord = RR'*cpts;
                        
                        physcoord(1,gp,elementCounter)     = coord(1);
                        physcoord(2,gp,elementCounter)     = coord(2);
                        physcoord(3,gp,elementCounter)     = coord(3);
                        
                        % Solve for first derivatives in global coordinates
                        dRdx = dxdxi\dR;
                        
                        cdRdx = dRdx(1,:);
                        cdRdy = dRdx(2,:);
                        cdRdz = dRdx(3,:);
                                                
                        % Calculate displacement values
                        dispcoord(1,gp,elementCounter) = RR'*tdisp(3*sctrx-2);
                        dispcoord(2,gp,elementCounter) = RR'*tdisp(3*sctrx-1);
                        dispcoord(3,gp,elementCounter) = RR'*tdisp(3*sctrx);

                        
                        strain11 = cdRdx*tdisp(3*sctrx-2);
                        strain22 = cdRdy*tdisp(3*sctrx-1);
                        strain33 = cdRdz*tdisp(3*sctrx);
                        strain12 = cdRdy*tdisp(3*sctrx-2) + cdRdx*tdisp(3*sctrx-1);
                        strain23 = cdRdy*tdisp(3*sctrx) + cdRdz*tdisp(3*sctrx-1);
                        strain31 = cdRdz*tdisp(3*sctrx-2) + cdRdx*tdisp(3*sctrx);                  
                           
                        epsiloncoord(:,gp,elementCounter) = [strain11;strain22;strain33;strain12;strain23;strain31];
                        wcoord(1,gp,elementCounter)= RR'*tdisp(3*sizeBasis+sctrx);   
                        gp = gp +1 ;
                    end
                end
            end
            elementCounter = elementCounter + 1;
        end
    end
end
msh_to_vtu_3dVM(physcoord,dispcoord,epsiloncoord,wcoord,[nx ny nz],vtuFile);
end

