function plotDispPhase2D_bimaterial(PHTelem,tdisp,sizeBasis,numberElements,geometry,controlPts,Mater,gaussCord)
% Evaluats the displacements and phase field values at the element vertices
% Displaying the displacements
% uses the PHT element structure
% for multipatch geometries
dim = geometry.dim;
p = geometry.p;
q = geometry.q;
vIEN = zeros(numberElements,4);
physcoord = zeros(4*numberElements,2);
dispcoord = zeros(4*numberElements,2);
wcoord = zeros(4*numberElements,1);
sigmacoord = zeros(4*numberElements,1);

fudge = 0;
nx = 2;
ny = 2;
px = linspace(-1+fudge,1-fudge, nx);
py = linspace(-1+fudge,1-fudge, ny);

%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u,dB_u] = bernstein_basis(px,p);
[B_v,dB_v] = bernstein_basis(py,q);

dBdu = zeros(nx, ny, (p+1)*(q+1));
dBdv = zeros(nx, ny, (p+1)*(q+1));
R = zeros(nx, ny, (p+1)*(q+1));

% The derivatives of the 2D Bernstein polynomials at Gauss points on the master element
basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        dBdu(:,:,basisCounter) = dB_u(:,i)*B_v(:,j)';
        dBdv(:,:,basisCounter) = B_u(:,i)*dB_v(:,j)';
        R(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
    end
end


elementCounter = 0;
for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            elementCounter =  elementCounter+1;
            
            vIEN(elementCounter,:) = [(elementCounter-1)*4+1:(elementCounter-1)*4+4];
            ximin = PHTelem{indexPatch}(i).vertex(1);
            ximax = PHTelem{indexPatch}(i).vertex(3);
            etamin = PHTelem{indexPatch}(i).vertex(2);
            etamax = PHTelem{indexPatch}(i).vertex(4);
            
            coordt = cell(ny,nx);
            dispmatx = zeros(ny,nx);
            dispmaty = zeros(ny,nx);
            stressvect = cell(ny,nx);
            stressVM = cell(ny,nx);
            wmatx = zeros(ny, nx);
            
            nument = size(PHTelem{indexPatch}(i).C,1);
            nodes = PHTelem{indexPatch}(i).nodes(1:nument);
            sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            dsctrx = reshape([2*sctrx-1;2*sctrx],1,2*nument);
            cpts = controlPts{indexPatch}(nodes,1:2);
            wgts = controlPts{indexPatch}(nodes,3);
            if (gaussCord{elementCounter}(end,2)<=20.0)
                C = Mater.C1;
            else
                C = Mater.C2;
            end
            
            for ii=1:nx
                for jj=1:ny
                    dRdx = (PHTelem{indexPatch}(i).C)*squeeze(dBdu(ii,jj,:))*2/(ximax-ximin);
                    dRdy = (PHTelem{indexPatch}(i).C)*squeeze(dBdv(ii,jj,:))*2/(etamax-etamin);
                    RR = (PHTelem{indexPatch}(i).C)*squeeze(R(ii,jj,:));
                    
                    RR = RR .* wgts;
                    dRdx = dRdx .* wgts;
                    dRdy = dRdy .* wgts;
                    
                    w_sum = sum(RR);
                    dw_xi = sum(dRdx);
                    dw_eta = sum(dRdy);
                    
                    dRdx = dRdx/w_sum - RR*dw_xi/w_sum^2;
                    dRdy = dRdy/w_sum - RR*dw_eta/w_sum^2;
                    RR = RR/w_sum;
                    
                    % multiply by the jacobian of the transformation from reference
                    % space to the parameter space
                    coord = RR'*cpts;
                    dR  = [dRdx';dRdy'];
                    dxdxi = dR*cpts;
                    
                    
                    % Solve for first derivatives in global coordinates
                    dRdx = dxdxi\dR;
                    
                    Bu = zeros(3,dim*nument);
                    for inode=1:nument
                        % Calculation of Bu
                        Bu(1,2*inode-1)=dRdx(1,inode);
                        Bu(2,2*inode)=dRdx(2,inode);
                        Bu(3,2*inode-1)=dRdx(2,inode);
                        Bu(3,2*inode)=dRdx(1,inode);
                    end
                    
                    coordt{jj,ii} = coord;
                    
                    
                    % Calculate tdisp values
                    dispmatx(jj,ii) = dispmatx(jj,ii) + RR'*tdisp(2*sctrx-1);
                    dispmaty(jj,ii) = dispmaty(jj,ii) + RR'*tdisp(2*sctrx);
                    
                    % Calculate stress values
                    stressvect{jj,ii} = C*Bu*tdisp(dsctrx);
                    stressVM{jj,ii} = sqrt(stressvect{jj,ii}(1)^2 - stressvect{jj,ii}(1)*stressvect{jj,ii}(2) + stressvect{jj,ii}(2)^2 +3*stressvect{jj,ii}(3)^2);
                    
                    %calculate phase field values
                    wmatx(jj,ii) = wmatx(jj,ii) + RR'*tdisp(2*sizeBasis+sctrx);
                    
                end
            end
            physcoord((elementCounter-1)*4+1,:) = coordt{1,1};
            physcoord((elementCounter-1)*4+2,:) = coordt{1,2};
            physcoord((elementCounter-1)*4+3,:) = coordt{2,2};
            physcoord((elementCounter-1)*4+4,:) = coordt{2,1};
            
            dispcoord((elementCounter-1)*4+1,:) = [dispmatx(1,1) dispmaty(1,1)];
            dispcoord((elementCounter-1)*4+2,:) = [dispmatx(1,2) dispmaty(1,2)];
            dispcoord((elementCounter-1)*4+3,:) = [dispmatx(2,2) dispmaty(2,2)];
            dispcoord((elementCounter-1)*4+4,:) = [dispmatx(2,1) dispmaty(2,1)];
            
            
            sigmacoord((elementCounter-1)*4+1, 1) = stressVM{1,1};
            sigmacoord((elementCounter-1)*4+2, 1) = stressVM{1,2}';
            sigmacoord((elementCounter-1)*4+3, 1) = stressVM{2,2}';
            sigmacoord((elementCounter-1)*4+4, 1) = stressVM{2,1}';
            
            wcoord((elementCounter-1)*4+1) = wmatx(1,1);
            wcoord((elementCounter-1)*4+2) = wmatx(1,2);
            wcoord((elementCounter-1)*4+3) = wmatx(2,2);
            wcoord((elementCounter-1)*4+4) = wmatx(2,1);
            
        end
    end
end

% Magnifcation factor for the tdisp plot
factor = 10;

plot3 = subplot(2,2,3);
cla(plot3)
trisurf(vIEN,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), log10(sigmacoord), 'EdgeColor','none','facecolor','interp')
view(0,90)
title('Displacements and \sigma_{VM}')
colorbar('vert')
colormap('parula')
hold on
axis tight
axis equal

plot4 = subplot(2,2,4);
cla(plot4)
factor = 0;
trisurf(vIEN,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), wcoord, 'facecolor','interp', 'EdgeColor','none')
view(0,90)
title('Phase Field')
colorbar('vert')
colormap('parula')
hold on
axis tight
axis equal
drawnow

numVertices=4*size(vIEN,1);
end
