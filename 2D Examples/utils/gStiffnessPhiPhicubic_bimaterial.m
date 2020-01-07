function [stiffkPhiPhi,Phi] = gStiffnessPhiPhicubic_bimaterial(PHTelem,sizeBasis,numberElements,dgdx,shape,Fract,Mater,volume,geometry,fenerg,gaussCord,tdisp)
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
            sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            ephi = tdisp(sctrx+2*sizeBasis)';
            indexCounter = indexCounter + nument^2;
            localkPhiPhi = zeros(nument,nument);
            localfPhi = zeros(nument,1);
            kgauss = 0;
            if (gaussCord{elementCounter}(end,2)<=20.0)
                cenerg = Fract.cenerg1;
                C = Mater.C1;
            else
                cenerg = Fract.cenerg2;
                C = Mater.C2;
            end
            
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss+1;
                    phigp = ephi*shape{elementCounter}(kgauss,:)';
                    ddphi = 6.0*phigp-Fract.s*(3*phigp - 1);
                    [~,Bphi,~] = strainGrad(dgdx{elementCounter},nument,nstress,dim,kgauss,C);
                    
                    senerg = fenerg(elementCounter, kgauss);
                    
                    % Calculation of kPhiPhi
                    localkPhiPhi = localkPhiPhi + cenerg*Fract.constl* (Bphi'*Bphi).*volume(elementCounter,kgauss);
                    localkPhiPhi = localkPhiPhi + (((cenerg/Fract.constl) + ddphi*senerg))...
                        .*(shape{elementCounter}(kgauss,:)'*shape{elementCounter}(kgauss,:)).*volume(elementCounter,kgauss);
                    localfPhi = localfPhi + ddphi*shape{elementCounter}(kgauss,:)'*senerg*volume(elementCounter,kgauss);
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
Phi = zeros(sizeBasis,1);
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
stiffkPhiPhi = sparse(II,JJ,S,sizeBasis,sizeBasis);
end

