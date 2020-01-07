function [Basis] = calcBasis1D4th(PHTelem,controlPts,Mesh)

[Gpts,GWts]=GaussQuad(Mesh.nGauss);
[B_u,dB_u,ddB_u] = bernstein_basis2nd(Gpts,Mesh.p);
numPatches = length(PHTelem);

% Initialize the volume, dgdx and shape arrays
volume = cell(1,numPatches);
dgdx = cell(1,numPatches);
d2gdx2 = cell(1,numPatches);
shape = cell(1,numPatches);
gaussCord = cell(1,numPatches);

for indexPatch = 1:length(PHTelem)
    
    volume(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    dgdx(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    d2gdx2(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    shape(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    gaussCord(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            umin = PHTelem{indexPatch}(i).vertex(1);
            umax = PHTelem{indexPatch}(i).vertex(2);
            scalefac = (umax - umin)/2;
            nument = size(PHTelem{indexPatch}(i).C,1);
            nodes = PHTelem{indexPatch}(i).nodes(1:nument);
            cpts = controlPts{indexPatch}(nodes,1);
            wgts = controlPts{indexPatch}(nodes,2);
            kgauss = 0;
            for igaussX=1:Mesh.nGauss
                kgauss = kgauss + 1;
                RR = (PHTelem{indexPatch}(i).C)*B_u(igaussX,:)';
                dRdx = (PHTelem{indexPatch}(i).C)*dB_u(igaussX,:)'*2/(umax-umin);
                ddRdx = (PHTelem{indexPatch}(i).C)*squeeze(ddB_u(igaussX,:))'*4/((umax-umin).^2);
                
                RR = RR .* wgts;
                dRdx = dRdx .* wgts;
                ddRdx = ddRdx.*wgts;
                
                w_sum = sum(RR);
                dw_xi = sum(dRdx);
                d2w_xi = sum(ddRdx);
                
                ddRdx = ddRdx/w_sum - (2*dRdx*dw_xi + RR*d2w_xi)/w_sum^2 + 2*RR*dw_xi^2/w_sum^3; %dxidxi derivative  
                dRdx = dRdx/w_sum - RR*dw_xi/w_sum^2;
                RR = RR/w_sum;
                
                dR  = dRdx';
                ddR = ddRdx';
                
                dxdxi = dR*cpts;
                d2xdxi2 = ddR*cpts;
                dxdxi2 = dxdxi(1,1)^2 ;
                   
                coord = RR'*cpts;
                Basis.gaussCord{indexPatch}{i}(kgauss,:) = coord(:);
                % Solve for first derivatives in global coordinates
                dR = dxdxi\dR;
                ddR = (ddR'-dR'*d2xdxi2')/(dxdxi2);
                
                J = det(dxdxi);
                ddR =ddR';
                Basis.volume{indexPatch}{i}(kgauss) = J*scalefac*GWts(igaussX);
                Basis.shape{indexPatch}{i}(kgauss,:) = RR(:);
                Basis.dgdx{indexPatch}{i}(kgauss,:) = dR(1,:);
                Basis.d2gdx2{indexPatch}{i}(kgauss,:) = ddR(1,:);
            end
        end
    end
end
end