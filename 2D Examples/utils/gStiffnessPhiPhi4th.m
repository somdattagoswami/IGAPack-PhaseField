function [stiffkPhiPhi,Phi] = gStiffnessPhiPhi4th(PHTelem,dimBasis,numberElements,dgdx,d2gdx2,shape,Fract,Mater,volume,fenerg,geometry)
% Assembles the stiffness matrix and RHS (Galerkin method)

dim = geometry.dim;
nstress = geometry.nstress;
ngaussX = geometry.ngaussX;
ngaussY = geometry.ngaussY;
kPhiPhi =cell(numberElements,1);
fPhi =cell(numberElements,1);
elementCounter = 0;
indexCounter = 0;
for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            elementCounter =  elementCounter+1;
            nument = size(PHTelem{indexPatch}(i).C,1);
            indexCounter = indexCounter + nument^2;
            localkPhiPhi = zeros(nument,nument);
            localfPhi = zeros(nument,1);
            kgauss = 0;
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss+1;  
                    [~,Bphi,~,Hu] = strainGradHu(dgdx{elementCounter},d2gdx2{elementCounter},nument,nstress,dim,kgauss,Mater.C);
                    senerg = fenerg(elementCounter, kgauss);
                    
                    % Calculation of kPhiPhi                    
                    localkPhiPhi = localkPhiPhi + 0.5*Fract.cenerg*Fract.constl* (Bphi'*Bphi).*volume(elementCounter,kgauss);
                    localkPhiPhi = localkPhiPhi + (((Fract.cenerg/Fract.constl) + 2.0*senerg))...
                        .*(shape{elementCounter}(kgauss,:)'*shape{elementCounter}(kgauss,:)).*volume(elementCounter,kgauss);  
                    localkPhiPhi = localkPhiPhi + (1/16)*Fract.cenerg*(Fract.constl)^3*(Hu'*Hu).*volume(elementCounter,kgauss);
                    localfPhi = localfPhi +2*shape{elementCounter}(kgauss,:)'*senerg*volume(elementCounter,kgauss);
                end %jgaus
            end %igaus
            fPhi{elementCounter} = localfPhi;
            kPhiPhi{elementCounter} = localkPhiPhi;
        end
    end
end
II = zeros(1,indexCounter);
JJ = zeros(1,indexCounter);
S = zeros(1, indexCounter);
Phi = zeros(dimBasis,1);
% Assembling the Stiffness Matrix
indexCounter = 0;
elementCounter = 0;
for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            elementCounter = elementCounter+1;
            nument = size(PHTelem{indexPatch}(i).C,1);
            sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);            

            
            II(indexCounter+1:indexCounter + nument^2) = repmat(sctrx,1,nument);
            JJ(indexCounter+1:indexCounter + nument^2) = reshape(repmat(sctrx,nument,1),1,nument^2);
            S(indexCounter+1:indexCounter + nument^2) = reshape(kPhiPhi{elementCounter},1,nument^2);
            Phi(sctrx) = Phi(sctrx) + fPhi{elementCounter};
            indexCounter = indexCounter + nument^2;
        end
    end
end
stiffkPhiPhi = sparse(II,JJ,S,dimBasis,dimBasis);
end

