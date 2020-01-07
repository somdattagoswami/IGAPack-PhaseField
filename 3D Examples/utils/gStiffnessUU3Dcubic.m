function [stiff,elemRef] = gStiffnessUU3Dcubic(numberElements,sctrxElem,dimBasis,dgdxCell,shapeCell,volumeCell,elemIndexCell,Fract,Mater,geometry,solPhi)
% Assembles the stiffness matrix and RHS (Galerkin method)
% Gauss points

dim = geometry.dim;
iCounter = 0;

%calculate the length of the the triplet array
indexSize = zeros(1,numberElements);
for iElem = 1:numberElements
    sctrx = sctrxElem{iElem};
    nument = length(sctrx);
    indexSize(iElem) = dim^2*nument^2;
end

%initialize the triplet arrays
lengthTriplet = sum(indexSize);
II = zeros(1,lengthTriplet);
JJ = zeros(1,lengthTriplet);
S = zeros(1,lengthTriplet);
elemRef = zeros(1,numberElements);
for iElem = 1:numberElements
    
        
    sctrx = sctrxElem{iElem};
    nument = length(sctrx);
    dsctrx = reshape([3*sctrx-2;3*sctrx-1;3*sctrx],1,dim*nument);
    ephi = solPhi(sctrx);
    localkUU = zeros(dim*nument,dim*nument);
    
    dgdx = dgdxCell{iElem};
    shape = shapeCell{iElem};
    volume = volumeCell{iElem};
    elemIndex = elemIndexCell{iElem};
    for ii=1:geometry.ngaussX
        for jj=1:geometry.ngaussY
            for kk = 1:geometry.ngaussZ
                Bu = zeros(geometry.nstress,geometry.dim*nument);
                %kgauss = kgauss+1;
                kgauss = kk+geometry.ngaussY*(jj-1)+geometry.ngaussX*geometry.ngaussY*(ii-1);
                phigp = shape(:,kgauss)'*ephi;
                if (phigp > geometry.threshPhi) && (elemIndex(1,2) < geometry.maxRefLevel)
                        elemRef(iElem)= elemIndex(1,1); % Obtaining the element numbers
                end
                
                Bu(1,1:3:3*nument-2) = dgdx(1,:,kgauss);
                Bu(2,2:3:3*nument-1) = dgdx(2,:,kgauss);
                Bu(3,3:3:3*nument) =  dgdx(3,:,kgauss);
                
                Bu(4,1:3:3*nument-2) = dgdx(2,:,kgauss);
                Bu(4,2:3:3*nument-1) = dgdx(1,:,kgauss);
                
                Bu(5,2:3:3*nument-1) = dgdx(3,:,kgauss);
                Bu(5,3:3:3*nument) = dgdx(2,:,kgauss);
                
                Bu(6,1:3:3*nument-2) = dgdx(3,:,kgauss);
                Bu(6,3:3:3*nument) = dgdx(1,:,kgauss);
                
                gphi = Fract.s*((1.0-phigp).^3 -(1.0-phigp).^2) + 3*(1.0-phigp).^2 - 2*(1.0-phigp).^3;
                
                % Calculation of kUU
                localkUU = localkUU+ gphi.*(Bu'*Mater.C*Bu).*volume(kgauss);
            end
        end %jgaus
    end %igaus
    II(iCounter+1:iCounter+indexSize(iElem)) = repmat(dsctrx,1,dim*nument);
    JJ(iCounter+1:iCounter+indexSize(iElem)) = reshape(repmat(dsctrx,dim*nument,1),1,dim^2*nument^2);
    S(iCounter+1:iCounter+indexSize(iElem)) = reshape(localkUU,1,dim^2*nument^2);
    iCounter = iCounter + indexSize(iElem);
    
end

stiff = sparse(II,JJ,S,dim*dimBasis,dim*dimBasis);


end

