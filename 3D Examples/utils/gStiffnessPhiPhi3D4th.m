function [stiffPhiPhi,RHSPhi,fenerg] = gStiffnessPhiPhi3D4th(numberElements,sctrxElem,sizeBasis,dgdxCell,d2gdx2Cell,shapeCell,volumeCell,Fract,Mater,geometry,fenerg,tdisp)
% Assembles the Phi part of the stiffness matrix and RHS (Galerkin method)
% Gauss points
RHSPhi = zeros(sizeBasis,1);
%calculate the length of the the triplet array
indexSize = zeros(1,numberElements);
for iElem = 1:numberElements
    sctrx = sctrxElem{iElem};
    nument = length(sctrx);
    indexSize(iElem) = nument^2;
end

%initialize the triplet arrays
lengthTriplet = sum(indexSize);
II = zeros(1,lengthTriplet);
JJ = zeros(1,lengthTriplet);
S = zeros(1,lengthTriplet);

iCounter = 0;

for iElem = 1:numberElements
    sctrx = sctrxElem{iElem};
    nument = length(sctrx);
    dsctrx = reshape([3*sctrx-2;3*sctrx-1;3*sctrx],1,geometry.dim*nument);
    dispElmt = tdisp(dsctrx);
    
    localkPhiPhi = zeros(nument,nument);
    localPhi = zeros(nument,1);
    kgauss = 0;
    dgdx = dgdxCell{iElem};
    d2gdx2 = d2gdx2Cell{iElem};
    volume = volumeCell{iElem};
    shape = shapeCell{iElem};
    BuPhi = zeros(geometry.dim,geometry.dim*nument);
    
    for ii=1:geometry.ngaussX
        for jj=1:geometry.ngaussY
            for kk = 1:geometry.ngaussZ
                kgauss = kgauss+1;
                % Calculation of Bphi
                BuPhi(1,1:3:3*nument-2) = dgdx(1,:,kgauss);
                BuPhi(2,2:3:3*nument-1) = dgdx(2,:,kgauss);
                BuPhi(3,3:3:3*nument) =  dgdx(3,:,kgauss);
                
                Hu(1,1:nument)= d2gdx2(1,:,kgauss);
                Hu(4,1:nument)= d2gdx2(4,:,kgauss);
                Hu(7,1:nument)= d2gdx2(6,:,kgauss);
                Hu(11,1:nument)= d2gdx2(4,:,kgauss);
                Hu(14,1:nument)= d2gdx2(2,:,kgauss);
                Hu(17,1:nument)= d2gdx2(5,:,kgauss);
                Hu(21,1:nument)= d2gdx2(6,:,kgauss);
                Hu(24,1:nument)= d2gdx2(5,:,kgauss);
                Hu(27,1:nument)= d2gdx2(3,:,kgauss);
                Hu(28,1:nument)= d2gdx2(4,:,kgauss);
                Hu(29,1:nument)= d2gdx2(1,:,kgauss);
                Hu(31,1:nument)= d2gdx2(2,:,kgauss);
                Hu(32,1:nument)= d2gdx2(4,:,kgauss);
                Hu(34,1:nument)= d2gdx2(5,:,kgauss);
                Hu(35,1:nument)= d2gdx2(6,:,kgauss);
                Hu(37,1:nument)= d2gdx2(6,:,kgauss);
                Hu(39,1:nument)= d2gdx2(1,:,kgauss);
                Hu(40,1:nument)= d2gdx2(5,:,kgauss);
                Hu(42,1:nument)= d2gdx2(4,:,kgauss);
                Hu(43,1:nument)= d2gdx2(3,:,kgauss);
                Hu(45,1:nument)= d2gdx2(6,:,kgauss);
                Hu(47,1:nument)= d2gdx2(6,:,kgauss);
                Hu(48,1:nument)= d2gdx2(4,:,kgauss);
                Hu(50,1:nument)= d2gdx2(5,:,kgauss);
                Hu(51,1:nument)= d2gdx2(2,:,kgauss);
                Hu(53,1:nument)= d2gdx2(3,:,kgauss);
                Hu(54,1:nument)= d2gdx2(5,:,kgauss); 
                
                strainElmt = BuPhi*dispElmt;
                strainPos = (sum(strainElmt) + abs(sum(strainElmt)))/2;
                senerg = 0.5*Mater.lamda*(strainPos)^2 + Mater.mu*sum(strainPos.^2);
                senerg = max(senerg,fenerg{iElem}(kgauss));
                fenerg{iElem}(kgauss) = senerg;
                
                
                % Calculation of kPhiPhi
                localkPhiPhi = localkPhiPhi + 0.5*Fract.cenerg*Fract.constl* (dgdx(:,:,kgauss)'*dgdx(:,:,kgauss)).*volume(kgauss);
                localkPhiPhi = localkPhiPhi + (((Fract.cenerg/Fract.constl) + 2.0*senerg))...
                    .*(shape(:,kgauss)*shape(:,kgauss)').*volume(kgauss);
                localkPhiPhi = localkPhiPhi + (1/16)*Fract.cenerg*(Fract.constl)^3*(Hu'*Hu).*volume(kgauss);
                localPhi = localPhi + 2 *shape(:,kgauss)*senerg*volume(kgauss);
            end % kguass
        end %jgauss
    end %igauss
    
    RHSPhi(sctrx) = RHSPhi(sctrx) + localPhi;
    
    II(iCounter+1:iCounter+indexSize(iElem)) = repmat(sctrx,1,nument);
    JJ(iCounter+1:iCounter+indexSize(iElem)) = reshape(repmat(sctrx,nument,1),1,nument^2);
    S(iCounter+1:iCounter+indexSize(iElem)) = reshape(localkPhiPhi,1,nument^2);
    iCounter = iCounter + indexSize(iElem);
end
stiffPhiPhi = sparse(II,JJ,S,sizeBasis,sizeBasis);

