function [PHTelem,controlPts,dimBasis] = initMeshTensile(geometry)
% Uniform refinement for cubic PHT splines

% Initial mesh
numPatches = geometry.numPatches;
L = geometry.L;
W = geometry.W;
p = geometry.p;
q = geometry.q;
numberElemU = geometry.numElemU;
numberElemV = geometry.numElemV;

% Initialize the PHT geometry on coarse mesh
controlPts = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
dimBasis = zeros(1, numPatches);
quadList = cell(numPatches,1);
% Divide the patches along the x direction
xVertices = linspace(0,L,numPatches+1);
for patchIndex = 1:numPatches
    
    % Set the dimensions of the patch
    patchMinX = xVertices(patchIndex);
    patchMaxX = xVertices(patchIndex+1);
    patchMinY = 0;
    patchMaxY = W;
    
    % Initialize geometry on coarsest mesh
    coefs(1:4,1,1) = [patchMinX;patchMinY;0;1];
    coefs(1:4,1,2) = [patchMinX;patchMaxY;0;1];
    coefs(1:4,2,1) = [patchMaxX;patchMinY;0;1];
    coefs(1:4,2,2) = [patchMaxX;patchMaxY;0;1];
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    
    nurbs = nrbmak(coefs,{knotU,knotV});
    p_init = nurbs.order(1)-1;
    q_init = nurbs.order(2)-1;
    
    % Refine into numberElemU by numberElemV knotspans
    knotU = linspace(0,1,numberElemU+1);
    knotV = linspace(0,1,numberElemV+1);
    
    numberElementsU = length(unique(knotU))-1;
    numberElementsV = length(unique(knotV))-1;
    
    
    % Increase polynomial order
    nurbs = nrbdegelev(nurbs,[p-p_init,q-q_init]);
    nurbs = nrbkntins(nurbs,{knotU(2:end-1) knotV(2:end-1)});
    % Repeat the knots to get C1 continuity
    nurbs = nrbkntins(nurbs,{knotU(2:end-1) knotV(2:end-1)});
    [controlPts{patchIndex},PHTelem{patchIndex},dimBasis(patchIndex),quadList{patchIndex}] = genControlPtsNoRep(nurbs,p,q,numberElementsU,numberElementsV);

end
