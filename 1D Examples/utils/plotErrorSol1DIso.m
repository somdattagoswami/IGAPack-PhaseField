function plotErrorSol1DIso(PHTelem,controlPts,sol0,p,RHSfun,Fract)
% Plot the errors in the solution
% Calculate the L2 norm and H1 seminorm of the error
numPlotPts = 101;
plotPts = linspace(-1,1,numPlotPts);

% Define the 1D Bernstein polynomials
[Bu,dBu] = bernstein_basis(plotPts, p);
numPatches = length(PHTelem);

for indexPatch = 1:numPatches
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            umin = PHTelem{indexPatch}(i).vertex(1);
            umax = PHTelem{indexPatch}(i).vertex(2);
            
            nument = size(PHTelem{indexPatch}(i).C,1); % Number of basis functions with support on current knotspan
            sctr = PHTelem{indexPatch}(i).nodes(1:nument);
            coord = zeros(1, numPlotPts);
            nodes = PHTelem{indexPatch}(i).nodes(1:nument);
            cpts = controlPts{indexPatch}(nodes,1);
            wgts = controlPts{indexPatch}(nodes,2);
            compSol = zeros(1, numPlotPts);
            compDeriv = zeros(1, numPlotPts);
            
            for ii=1:numPlotPts
                RR = (PHTelem{indexPatch}(i).C)*Bu(ii,:)';
                dRdx = (PHTelem{indexPatch}(i).C)*dBu(ii,:)'*2/(umax-umin);
                RR = RR .* wgts;
                dRdx = dRdx .* wgts;
                w_sum = sum(RR);
                dw_xi = sum(dRdx);
                dRdx = dRdx/w_sum - RR*dw_xi/w_sum^2;
                RR = RR/w_sum;
                dR  = dRdx';
                dxdxi = dR*cpts;
                coord(ii) = RR'*cpts;
                dR = dxdxi\dR;
                
                % Calculate displacement values
                compSol(ii) = RR'*sol0(sctr);
                
                % Calculate the error in stress values
                compDeriv(ii) = dR*sol0(sctr);
            end
            if nargin<6
                [analSol] = RHSfun(coord);
            else
                [analSol] = RHSfun(coord,Fract);
            end
            plot(coord,compSol,'-r','linewidth',3)
            hold on
            plot(coord,analSol,'-b','linewidth',3)
            hold on
            plot(coord, zeros(1,numPlotPts),'-k')
            hold on
        end
    end
end
title('Error in the solution')
drawnow
end

