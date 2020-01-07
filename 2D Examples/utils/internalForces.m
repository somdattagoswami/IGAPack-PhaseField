function fenerg = internalForces(PHTelem,dgdx,tdisp,geometry,Mater,fenerg)
% Compute the decomposed stress

dim = geometry.dim;
nstress = geometry.nstress;

elementCounter = 0;

for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            elementCounter =  elementCounter+1;
            nument = size(PHTelem{indexPatch}(i).C,1);
            sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            dsctrx = reshape([2*sctrx-1;2*sctrx],1,dim*nument);
            dispElmt = tdisp(dsctrx);

            % Calculate strains and stresses at integration points
            kgauss=0;
            for ii=1:geometry.ngaussX
                for jj=1:geometry.ngaussY
                    kgauss=kgauss+1;
                                            
                    [Bu,~,~]=strainGrad(dgdx{elementCounter},nument,nstress,dim,kgauss,Mater.C);
                    strainElmt = Bu*dispElmt;
                    [positive_elast]=localDStress2(strainElmt(1),strainElmt(2),strainElmt(3)/2,Mater.lamda,Mater.mu);
    
                    senerg = max(positive_elast,fenerg(elementCounter,kgauss));
                    fenerg(elementCounter,kgauss) = senerg;
    
                end %igaus
            end %jgaus
        end
    end
end

end %endfunction