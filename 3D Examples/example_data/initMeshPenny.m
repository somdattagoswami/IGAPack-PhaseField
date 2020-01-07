function [PHTelem,controlPts,dimBasis,numberElements,octupleList] = initMeshPenny(geometry)
% Non-uniform refinement for cubic PHT splines

% Initial mesh
numPatches = geometry.numPatches;
L = geometry.L;
W = geometry.W;
H = geometry.H;
numberElemU = geometry.numElemU;
numberElemV = geometry.numElemV;
numberElemW = geometry.numElemW;

% Initialize the PHT geometry on coarse mesh
controlPts = cell(numPatches,1);
PHTelem = cell(numPatches,1);
dimBasis = zeros(1,numPatches);
octupleList = cell(numPatches,1);

% Divide the patches along the x direction
xVertices = linspace(-L/2,L/2,numPatches+1);
numberElements = 0;
for patchIndex = 1:numPatches
    
    % Set the dimensions of the patch
    patchMinX = xVertices(patchIndex);
    patchMaxX = xVertices(patchIndex+1);
    patchMinY = -W/2;
    patchMaxY = W/2;
    patchMinZ = -H/2;
    patchMaxZ = H/2;
    
    % Initialize geometry on coarsest mesh
    coefs(1:3,1,1,1) = [patchMinX; patchMinY; patchMinZ];
    coefs(1:3,1,2,1) = [patchMinX; patchMaxY; patchMinZ];
    coefs(1:3,2,1,1) = [patchMaxX; patchMinY; patchMinZ];
    coefs(1:3,2,2,1) = [patchMaxX; patchMaxY; patchMinZ];
    coefs(1:3,1,1,2) = [patchMinX; patchMinY; patchMaxZ];
    coefs(1:3,1,2,2) = [patchMinX; patchMaxY; patchMaxZ];
    coefs(1:3,2,1,2) = [patchMaxX; patchMinY; patchMaxZ];
    coefs(1:3,2,2,2) = [patchMaxX; patchMaxY; patchMaxZ];
    coefs(4,1,1,1) = 1;
    coefs(4,1,2,1) = 1;
    coefs(4,2,1,1) = 1;
    coefs(4,2,2,1) = 1;
    coefs(4,1,1,2) = 1;
    coefs(4,1,2,2) = 1;
    coefs(4,2,1,2) = 1;
    coefs(4,2,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    knotW = [0 0 1 1];
    nurbs = nrbmak(coefs,{knotU,knotV,knotW});
    
    % Refine into numberElemU by numberElemV knotspans
    knotU = linspace(0,1,numberElemU+1);
    knotV = linspace(0,1,numberElemV+1);
    knotW = linspace(0,1,numberElemW+1);
    
    geometry.numElmtU = length(unique(knotU))-1;
    geometry.numElmtV = length(unique(knotV))-1;
    geometry.numElmtW = length(unique(knotW))-1;
    
    % Increase polynomial order
    p = geometry.p;
    q = geometry.q;
    r = geometry.r;
    nurbs = nrbdegelev(nurbs,[p,q,r]-(nurbs.order-1));
    nurbs = nrbkntins(nurbs,{knotU(2:end-1) knotV(2:end-1) knotW(2:end-1)});
    
    % Repeat the knots to get C1 continuity
    nurbs = nrbkntins(nurbs,{knotU(2:end-1) knotV(2:end-1) knotW(2:end-1)});
    
    [PHTelem{patchIndex},controlPts{patchIndex},dimBasis(patchIndex),octupleList{patchIndex}] = genControlPts3D(nurbs,geometry);
    numberElements = numberElements + 8*size(octupleList{patchIndex}(:,1),1);
end
end
