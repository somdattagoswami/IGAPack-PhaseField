function [stiffPhiPhi,RHSPhi,fenerg] = gStiffnessPhiPhi3Dcubic(numberElements,sctrxElem,sizeBasis,dgdxCell,shapeCell,volumeCell,Fract,Mater,geometry,fenerg,solU,solPhi)
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
    dispElmt = solU(dsctrx);
    ephi = solPhi(sctrx)';
    
    localkPhiPhi = zeros(nument,nument);
    localPhi = zeros(nument,1);
    kgauss = 0;
    dgdx = dgdxCell{iElem};    
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
                
                strainElmt = BuPhi*dispElmt;
                strainPos = (sum(strainElmt) + abs(sum(strainElmt)))/2;
                senerg = 0.5*Mater.lamda*(strainPos)^2 + Mater.mu*sum(strainPos.^2);
                senerg = max(senerg,fenerg{iElem}(kgauss));
                fenerg{iElem}(kgauss) = senerg;
                phigp = ephi*shape(:,kgauss);
                ddphi = 6.0*phigp - Fract.s*(3*phigp - 1);
                
                % Calculation of kPhiPhi
                localkPhiPhi = localkPhiPhi + Fract.cenerg*Fract.constl* (dgdx(:,:,kgauss)'*dgdx(:,:,kgauss)).*volume(kgauss);
                localkPhiPhi = localkPhiPhi + (((Fract.cenerg/Fract.constl) + ddphi*senerg))...
                    .*(shape(:,kgauss)*shape(:,kgauss)').*volume(kgauss);
                localPhi = localPhi + ddphi*shape(:,kgauss)*senerg*volume(kgauss);
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

