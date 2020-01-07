function [stiffkUU,rhs,elemRef] = gStiffnessUU1DMark(PHTelem,Basis,Fract,Mater,tdisp,Mesh,degFun,RHSFun)

dim = Mesh.dim;
numPatches = length(PHTelem);
rhs = zeros(dim*Mesh.sizeBasis,1);
elemRef = cell(numPatches,1);
indexCounter = 0;
for indexPatch = 1:length(PHTelem)
    elemRef(indexPatch) = {zeros(1,length(PHTelem{indexPatch}))};
    for i=1:length(PHTelem{indexPatch})        
        if isempty(PHTelem{indexPatch}(i).children)
            nument = size(PHTelem{indexPatch}(i).C,1);
            sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            ephi = tdisp(sctrx + Mesh.sizeBasis)';
            localkUU = zeros(dim*nument,dim*nument);
            localrhs = zeros(dim*nument,1); %local RHS vector
            kgauss = 0;
            for ii=1:Mesh.nGauss                
                    kgauss = kgauss+1;
                    phigp = ephi*Basis.shape{indexPatch}{i}(kgauss,:)';
                    if (phigp > Mesh.threshPhi) && (PHTelem{indexPatch}(i).level <= Mesh.maxRefLevel)
                        elemRef{indexPatch}(i) = 1;
                    end
                    
                    Bu = Basis.dgdx{indexPatch}{i}(kgauss,:); 
                    D = Mater.E*Bu;
                    coord = Basis.gaussCord{indexPatch}{i}(kgauss,:);
                    P = RHSFun(coord);
                    % Calculation of kUU
                    gphi = degFun(phigp);
                    localkUU = localkUU+(gphi+Fract.constk).*(Bu'*D).*Basis.volume{indexPatch}{i}(kgauss);
                    localrhs = localrhs + Basis.shape{indexPatch}{i}(kgauss,:)'*P*Basis.volume{indexPatch}{i}(kgauss);
            end %igaus
            II(indexCounter+1:indexCounter+dim^2*nument^2) = repmat(sctrx,1,dim*nument);
            JJ(indexCounter+1:indexCounter+dim^2*nument^2) = reshape(repmat(sctrx,dim*nument,1),1,dim^2*nument^2);
            S(indexCounter+1:indexCounter+dim^2*nument^2) = reshape(localkUU,1,dim^2*nument^2);
            indexCounter = indexCounter +dim^2*nument^2;
            rhs(sctrx) = rhs(sctrx) + localrhs;
        end
    end
end
stiffkUU = sparse(II,JJ,S,dim*Mesh.sizeBasis,dim*Mesh.sizeBasis);
end

