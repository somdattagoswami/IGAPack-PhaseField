function [stiffkPhiPhi,Phi] = gStiffnessPhiPhi(PHTelem,dimBasis,numberElements,dgdx,shape,Fract,Mater,volume,geometry,fenerg)
% Assembles the stiffness matrix

dim = geometry.dim;
nstress = geometry.nstress;

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
            for ii=1:geometry.ngaussX
                for jj=1:geometry.ngaussY
                    kgauss = kgauss+1;  
                    [~,Bphi,~]=strainGrad(dgdx{elementCounter},nument,nstress,dim,kgauss,Mater.C);
                    
                    senerg = fenerg(elementCounter, kgauss);

                  
                    % Calculation of kPhiPhi                    
                    localkPhiPhi = localkPhiPhi + Fract.cenerg*Fract.constl* (Bphi'*Bphi).*volume(elementCounter,kgauss);
                    localkPhiPhi = localkPhiPhi + (((Fract.cenerg/Fract.constl) + 2.0*senerg))...
                        .*(shape{elementCounter}(kgauss,:)'*shape{elementCounter}(kgauss,:)).*volume(elementCounter,kgauss);   
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

