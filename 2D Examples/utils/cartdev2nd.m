function[shape,dgdx,d2gdx2,volume,gaussCord,numberElements]=cartdev2nd(PHTelem,controlPts,geometry)
%compute the derivatives using the PHTelem structure
%Gauss points

p = geometry.p;
q = geometry.q;

ngaussX = geometry.ngaussX;
ngaussY = geometry.ngaussY;
mgauss = ngaussX*ngaussY;
dim = geometry.dim;
[GptsX,GWtsX]=GaussQuad(ngaussX);
[GptsY,GWtsY]=GaussQuad(ngaussY);


%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u,dB_u,ddB_u] = bernstein_basis2nd(GptsX,p);
[B_v,dB_v,ddB_v] = bernstein_basis2nd(GptsY,q);

ddBdu2 = zeros(ngaussX, ngaussY, (p+1)*(q+1));
ddBdv2 = zeros(ngaussX, ngaussY, (p+1)*(q+1));
ddBduv = zeros(ngaussX, ngaussY, (p+1)*(q+1));
dBdu = zeros(ngaussX, ngaussY, (p+1)*(q+1));
dBdv = zeros(ngaussX, ngaussY, (p+1)*(q+1));
R = zeros(ngaussX, ngaussY, (p+1)*(q+1));

% The derivatives of the 2D Bernstein polynomials at Gauss points on the master element
basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        R(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
        dBdu(:,:,basisCounter) = dB_u(:,i)*B_v(:,j)';
        dBdv(:,:,basisCounter) = B_u(:,i)*dB_v(:,j)';
        ddBdu2(:,:,basisCounter) = ddB_u(:,i)*B_v(:,j)';
        ddBdv2(:,:,basisCounter) = B_u(:,i)*ddB_v(:,j)';
        ddBduv(:,:,basisCounter) =dB_u(:,i)*dB_v(:,j)';
    end
end

elementCounter = 0;
for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            elementCounter = elementCounter + 1;
        end
    end
end
numberElements = elementCounter;

% Initialize the volume, dgdx and shape arrays
volume = zeros(numberElements,mgauss);
dgdx = cell(numberElements,1);
d2gdx2 = cell(numberElements,1);
shape = cell(numberElements,1);
gaussCord = cell(numberElements,1);
elementCounter = 0;
for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            elementCounter = elementCounter + 1;
            
            xmin = PHTelem{indexPatch}(i).vertex(1);
            xmax = PHTelem{indexPatch}(i).vertex(3);
            ymin = PHTelem{indexPatch}(i).vertex(2);
            ymax = PHTelem{indexPatch}(i).vertex(4);
            
            % The jacobian of the transformation from [-1,1]x[-1,1] to [xmin,xmax]x[ymin,ymax]
            scalefac = (xmax - xmin)*(ymax - ymin)/4;
            nument = size(PHTelem{indexPatch}(i).C,1);
            nodes = PHTelem{indexPatch}(i).nodes(1:nument);
            cpts = controlPts{indexPatch}(nodes,1:2);
            wgts = controlPts{indexPatch}(nodes,3);
            deriv = zeros(mgauss,dim,nument);
            dderiv = zeros(mgauss,dim+1,nument);
            funcVal = zeros(mgauss,nument);
            kgauss = 0;
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss + 1;
                    ddRdx2 = (PHTelem{indexPatch}(i).C)*squeeze(ddBdu2(ii,jj,:))*4/((xmax-xmin).^2);
                    ddRdy2 = (PHTelem{indexPatch}(i).C)*squeeze(ddBdv2(ii,jj,:))*4/((ymax-ymin).^2);
                    ddRdxdy = (PHTelem{indexPatch}(i).C)*squeeze(ddBduv(ii,jj,:))*4/((ymax-ymin)*(xmax-xmin));
                    dRdx = (PHTelem{indexPatch}(i).C)*squeeze(dBdu(ii,jj,:))*2/(xmax-xmin);
                    dRdy = (PHTelem{indexPatch}(i).C)*squeeze(dBdv(ii,jj,:))*2/(ymax-ymin);
                    RR = (PHTelem{indexPatch}(i).C)*squeeze(R(ii,jj,:));
                    
                    RR = RR .* wgts;
                    dRdx = dRdx .* wgts;
                    dRdy = dRdy .* wgts;
                    ddRdx2 = ddRdx2.*wgts;
                    ddRdy2 = ddRdy2.*wgts;
                    ddRdxdy = ddRdxdy.*wgts;
                    
                    w_sum = sum(RR);
                    dw_xi = sum(dRdx);
                    dw_eta = sum(dRdy);
                    d2w_xi = sum(ddRdx2);
                    d2w_eta = sum(ddRdy2);
                    d2w_xieta = sum(ddRdxdy);
                    
                    ddRdx2 = ddRdx2/w_sum - (2*dRdx*dw_xi + RR*d2w_xi)/w_sum^2 + 2*RR*dw_xi^2/w_sum^3; %dxidxi derivative
                    ddRdy2 = ddRdy2/w_sum - (2*dRdy*dw_eta + RR*d2w_eta)/w_sum^2 + 2*RR*dw_eta^2/w_sum^3; %detadeta derivative
                    ddRdxdy = ddRdxdy/w_sum - (dRdx*dw_eta + dRdy*dw_xi + RR*d2w_xieta)/w_sum^2 + 2*RR*dw_xi*dw_eta/w_sum^3; %dxideta derivative
                    dRdx = dRdx/w_sum - RR*dw_xi/w_sum^2;
                    dRdy = dRdy/w_sum - RR*dw_eta/w_sum^2;
                    RR = RR/w_sum;
                    
                    % multiply by the jacobian of the transformation from reference
                    % space to the parameter space
                    dR  = [dRdx';dRdy'];
                    ddR = [ddRdx2';ddRdy2';ddRdxdy'];
                    coord = RR'*cpts;
                    %plot(coord(1),coord(2),'.r')
                    %hold on
                    gaussCord{elementCounter}(kgauss,:) = coord(:);
                    dxdxi = dR*cpts;
                    d2xdxi2 = ddR*cpts; % Set up the second derivatives matrix and the matrix of squared first derivatives
                    dxdxi2 = [dxdxi(1,1)^2, dxdxi(1,2)^2, dxdxi(1,1)*dxdxi(1,2);...
                              dxdxi(2,1)^2, dxdxi(2,2)^2, dxdxi(2,1)*dxdxi(2,2), ;...
                              2*dxdxi(1,1)*dxdxi(2,1), dxdxi(1,2)*dxdxi(2,2),(dxdxi(1,1)*dxdxi(2,2))+(dxdxi(1,2)*dxdxi(2,1))] ;
                    % Solve for first derivatives in global coordinates
                    dR = dxdxi\dR;
                    ddR = (ddR'-dR'*d2xdxi2')/(dxdxi2);
                    J = det(dxdxi);
                    ddR =ddR';
                    volume(elementCounter,kgauss) = J*scalefac*GWtsX(ii).*GWtsY(jj);
                    for idime=1:dim
                        for inode=1:nument
                            funcVal(kgauss,inode) = RR(inode);
                            deriv(kgauss,idime,inode) = dR(idime,inode);
                            dderiv(kgauss,idime,inode) = ddR(idime,inode);
                        end
                    end
                    for inode=1:nument
                        dderiv(kgauss,dim+1,inode) = ddR(dim+1,inode);
                    end
                end
            end
            shape{elementCounter} = funcVal;
            dgdx{elementCounter} = deriv;
            d2gdx2{elementCounter} = dderiv;
        end
    end
end

