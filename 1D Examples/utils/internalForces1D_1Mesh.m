function [fenerg] = internalForces1D_1Mesh(PHTelem,dgdx,tdisp,Mater,fenerg,nGauss)

for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            nument = size(PHTelem{indexPatch}(i).C,1);
            sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            dispElmt = tdisp(sctrx);
            % Calculate strains and stresses at integration points
            kgauss=0;
            for ii=1:nGauss                
                    kgauss=kgauss+1;
                    Bu = dgdx{indexPatch}{i}(kgauss,:);
                    strainElmt = Bu*dispElmt;
                    stress =  Mater.E*strainElmt;
                    senerg = 0.5*stress* strainElmt;
                    senerg = max(senerg,fenerg{indexPatch}{i}(kgauss));
                    fenerg{indexPatch}{i}(kgauss) = senerg;               
            end %iigaus
        end
    end
end

