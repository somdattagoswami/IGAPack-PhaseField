function [PHTelem,controlPts,dirichlet,Mesh] = initMesh(Mesh)
% Non-uniform refinement for cubic PHT splines
% Initial mesh
numPatches = Mesh.numPatches;
L = Mesh.L;

% Initialize the PHT geometry on coarse mesh
controlPts = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
dimBasis = zeros(1, numPatches);
tupleList = cell(numPatches,1);
% Divide the patches along the x direction
xVertices = linspace(-L,L,numPatches+1);
numberElements = 0;
for patchIndex = 1:numPatches
    
    % Set the dimensions of the patch
    patchMinX = xVertices(patchIndex);
    patchMaxX = xVertices(patchIndex+1);
    
    % Initialize geometry on coarsest mesh
    coefs(1:4,1) = [patchMinX;0;0;1];
    coefs(1:4,2) = [patchMaxX;0;0;1];
    
    knotU = [0 0 1 1];
    nurbs = nrbmak(coefs,knotU);
    p_init = nurbs.order(1)-1;
    % Refine into numberElemU knotspans
    knotU = linspace(0,1,Mesh.numElem+1);
    numberElementsU = length(unique(knotU))-1;
    
    % Increase polynomial order
    nurbs = nrbdegelev(nurbs,Mesh.p-p_init);
    nurbs = nrbkntins(nurbs,knotU(2:end-1));
    % Repeat the knots to get C1 continuity
    nurbs = nrbkntins(nurbs,knotU(2:end-1));
    [controlPts{patchIndex},PHTelem{patchIndex},dimBasis(patchIndex), ...
        tupleList{patchIndex}] = genControlPts1D(nurbs,Mesh.p,numberElementsU);
    numberElements = numberElements + 2*size(tupleList{patchIndex},1);
    [dirichlet] = initialCrack1D(PHTelem{patchIndex});
end
Mesh.numberElements = numberElements;
Mesh.tupleList = tupleList;
Mesh.dimBasis = dimBasis;
