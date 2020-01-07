function fenerg = internalForces_bimaterial(PHTelem,dgdx,tdisp,geometry,Mater,fenerg,gaussCord)
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
            
            if (gaussCord{elementCounter}(end,2)<=20.0)
                C = Mater.C1;
                lamda = Mater.lamda1;
                mu = Mater.mu1;
            else
                C = Mater.C2;
                lamda = Mater.lamda2;
                mu = Mater.mu2;
            end
            
            % Calculate strains and stresses at integration points
            kgauss=0;
            
            for ii=1:geometry.ngaussX
                for jj=1:geometry.ngaussY
                    kgauss=kgauss+1;
                    
                    [Bu,~,~]=strainGrad(dgdx{elementCounter},nument,nstress,dim,kgauss,C);
                    strainElmt = Bu*dispElmt;
                    [positive_elast]=localDStress2(strainElmt(1),strainElmt(2),strainElmt(3)/2,lamda,mu);
                    
                    senerg = max(positive_elast,fenerg(elementCounter,kgauss));
                    fenerg(elementCounter,kgauss) = senerg;
                    
                end %igaus
            end %jgaus
        end
    end
end

end %endfunction