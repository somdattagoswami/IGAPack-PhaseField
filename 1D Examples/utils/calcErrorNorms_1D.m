function  [l2relerr]=calcErrorNorms_1D(sol0,PHTelem,controlPts,p,RHSfun)
% Calculates the errors in the solution
numGaussX = p+1;
[Gpts,GWts]=GaussQuad(numGaussX);
[B_u,dB_u] = bernstein_basis(Gpts,p);
numPatches = length(PHTelem);

l2norm = 0;
l2relerr = 0;

for patchIndex = 1:numPatches
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            
            umin = PHTelem{patchIndex}(i).vertex(1);
            umax = PHTelem{patchIndex}(i).vertex(2);
            
            nument = size(PHTelem{patchIndex}(i).C,1); %number of basis functions with support on current knotspan
            sctr = PHTelem{patchIndex}(i).nodesGlobal(1:nument);
            nodes = PHTelem{patchIndex}(i).nodes(1:nument);
            cpts = controlPts{patchIndex}(nodes,1);
            wgts = controlPts{patchIndex}(nodes,2);
            scalefac = (umax - umin)/2;
            
            for ii=1:numGaussX  
                R = (PHTelem{patchIndex}(i).C)*B_u(ii,:)';
                dRdx = (PHTelem{patchIndex}(i).C)*dB_u(ii,:)'*2/(umax-umin);
                R = R .* wgts;
                dRdx = dRdx .* wgts;                
                w_sum = sum(R);
                dw_xi = sum(dRdx);
                dRdx = dRdx/w_sum - R*dw_xi/w_sum^2;
                dR = dRdx';
                R = R/w_sum;
                dxdxi = dR*cpts;
                coord = R'*cpts;                
                J = det(dxdxi);
                % Calculate displacement values
                compSol = R'*sol0(sctr);
                % Analytical solution
                [analSol] = RHSfun(coord);     
                l2norm = l2norm + analSol^2*J*scalefac*GWts(ii);                
                l2relerr = l2relerr + (analSol-compSol)^2*J*scalefac*GWts(ii);                
            end            
        end
    end
end
l2relerr = sqrt(l2relerr)/sqrt(l2norm);
end
