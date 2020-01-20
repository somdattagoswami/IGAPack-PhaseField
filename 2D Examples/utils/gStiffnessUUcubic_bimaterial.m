function [stiffkUU, elemRef] = gStiffnessUUcubic_bimaterial(PHTelem,sizeBasis,numberElements,dgdx,shape,gaussCord,Fract,Mater,volume,tdisp,geometry)
% Assembles the stiffness matrix and RHS (Galerkin method)

dim = geometry.dim;
nstress = geometry.nstress;
numPatches = geometry.numPatches;
ngaussX = geometry.ngaussX;
ngaussY = geometry.ngaussY;
kUU =cell(numberElements,1);
elemRef = cell(numPatches,1);
elementCounter = 0;
indexCounter = 0;

for indexPatch = 1:length(PHTelem)
    
    elemRef{indexPatch} = zeros(1, length(PHTelem{indexPatch}));
    
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            elementCounter =  elementCounter+1;
            nument = size(PHTelem{indexPatch}(i).C,1);
            sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            ephi = tdisp(sctrx+2*sizeBasis)';
            indexCounter = indexCounter + dim^2*nument^2;
            localkUU = zeros(dim*nument,dim*nument);
            kgauss = 0;
            if (gaussCord{elementCounter}(end,2)<=20.0)
                C = Mater.C1;
            else
                C = Mater.C2;
            end
            
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss+1;
                    
                    phigp = ephi*shape{elementCounter}(kgauss,:)';
                    %mark element for refinement if phase field value > 0.5
                    if (phigp > geometry.threshPhi) && (PHTelem{indexPatch}(i).level < geometry.maxRefLevel)
                        elemRef{indexPatch}(i)=1;
                    end
                    
                    [Bu,~,D]=strainGrad(dgdx{elementCounter},nument,nstress,dim,kgauss,C);
                    gphi = Fract.s*((1.0-phigp).^3 -(1.0-phigp).^2) + 3*(1.0-phigp).^2 - 2*(1.0-phigp).^3;
                    
                    % Calculation of kUU
                    localkUU=localkUU + gphi.*(Bu'*D).*volume(elementCounter,kgauss);
                end %jgaus
            end %igaus
            kUU{elementCounter} = localkUU;
            %norm(localkUU)
        end
    end
end
II = zeros(1,indexCounter);
JJ = zeros(1,indexCounter);
S = zeros(1, indexCounter);

% Assembling the Stiffness Matrix
indexCounter = 0;
elementCounter = 0;
for indexPatch = 1:numPatches
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            elementCounter = elementCounter+1;
            nument = size(PHTelem{indexPatch}(i).C,1);
            sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            dsctrx = reshape([2*sctrx-1;2*sctrx],1,dim*nument);            
            II(indexCounter+1:indexCounter+dim^2*nument^2) = repmat(dsctrx,1,dim*nument);
            JJ(indexCounter+1:indexCounter+dim^2*nument^2) = reshape(repmat(dsctrx,dim*nument,1),1,dim^2*nument^2);
            S(indexCounter+1:indexCounter+dim^2*nument^2) = reshape(kUU{elementCounter},1,dim^2*nument^2);
            indexCounter = indexCounter +dim^2*nument^2;
        end
    end
end
stiffkUU = sparse(II,JJ,S,dim*sizeBasis,dim*sizeBasis);
end

