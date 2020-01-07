function [stiffkPhiPhi,Phi] = gStiffnessPhiPhi1D4th(PHTelem,Basis,Mesh,Fract,fenerg,tdisp,degFunDeriv)
ngaussX = Mesh.nGauss;
Phi = zeros(Mesh.sizeBasis,1);
indexCounter = 0;
for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)            
            nument = size(PHTelem{indexPatch}(i).C,1);
            sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            ephi = tdisp(sctrx+Mesh.sizeBasis)';
            localkPhiPhi = zeros(nument,nument);
            localfPhi = zeros(nument,1);
            kgauss = 0;
            for ii=1:ngaussX                
                    kgauss = kgauss+1;
                    Bphi = Basis.dgdx{indexPatch}{i}(kgauss,:);
                    Hu = Basis.d2gdx2{indexPatch}{i}(kgauss,:);
                    senerg = fenerg{indexPatch}{i}(kgauss);
                    phigp = ephi*Basis.shape{indexPatch}{i}(kgauss,:)';
                    ddphi = degFunDeriv(phigp);

                    % Calculation of kPhiPhi                    
                    localkPhiPhi = localkPhiPhi + Fract.cenerg*Fract.constl*...
                        (Bphi'*Bphi).*Basis.volume{indexPatch}{i}(kgauss);
                    localkPhiPhi = localkPhiPhi + (((Fract.cenerg/Fract.constl) + ddphi*senerg))...
                        .*(Basis.shape{indexPatch}{i}(kgauss,:)'*Basis.shape{indexPatch}{i}(kgauss,:))...
                        .*Basis.volume{indexPatch}{i}(kgauss); 
                    localkPhiPhi = localkPhiPhi + (1/16)*Fract.cenerg*(Fract.constl)^3*...
                        (Hu'*Hu).*Basis.volume{indexPatch}{i}(kgauss);
                    localfPhi = localfPhi + ddphi*Basis.shape{indexPatch}{i}(kgauss,:)'...
                        *senerg*Basis.volume{indexPatch}{i}(kgauss);                
            end %igaus
            II(indexCounter+1:indexCounter + nument^2) = repmat(sctrx,1,nument);
            JJ(indexCounter+1:indexCounter + nument^2) = reshape(repmat(sctrx,nument,1),1,nument^2);
            S(indexCounter+1:indexCounter + nument^2) = reshape(localkPhiPhi,1,nument^2);
            Phi(sctrx) = Phi(sctrx) + localfPhi;
            indexCounter = indexCounter + nument^2;
        end
    end
end
stiffkPhiPhi = sparse(II,JJ,S,Mesh.sizeBasis,Mesh.sizeBasis);
end

