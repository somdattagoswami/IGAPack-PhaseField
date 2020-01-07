function [controlPts,PHTelem,dimBasis,tupleList] = genControlPts1D(nrb,p,numElemU)

dim = 1;
toleq = 1e-10;
knotU = nrb.knots;
coefs = nrb.coefs;
basisU = length(knotU)-p-1;
dimBasis = basisU;
controlPts = zeros(dimBasis,dim+1);
numElements = numElemU;
tupleList = zeros(numElements,2^dim);

for i=1:basisU % for each node in the x direction
    controlPts(i,:) = [coefs(1,i)./coefs(4,i), coefs(4,i)];
end %for i

PHTelem = struct;
% Initialize the neighbor connectivity lists
PHTelem.neighbor_left = [];
PHTelem.neighbor_right = [];

%loop through each element and compute the element-node connectivities
elementCounter = 0;

for i=1:length(knotU)-1
    if (knotU(i+1)>knotU(i)+toleq)  %the knotspan has non-zero area
        elementCounter = elementCounter + 1;
        PHTelem(elementCounter).parent = [];
        PHTelem(elementCounter).children = [];
        PHTelem(elementCounter).vertex = [knotU(i), knotU(i+1)];
        
        %we add the nodes from i-p...i in the u direction
        currow = (i-p):i;
        
        PHTelem(elementCounter).nodes=currow;
        PHTelem(elementCounter).level = 0;
        tupleList(elementCounter,:) = numElements + (2*(elementCounter-1)+1:2*elementCounter);
        
    end
end
% loop through each element and compute the neighbor lists and Bezier
% extraction operators
knotU = nrb.knots;
[C_u, ~] = bezierExtraction(knotU,p);
indexMatrix = 1:numElements;
for i=1:numElemU
    elementIndex = indexMatrix(i);
    PHTelem(elementIndex).C = C_u(:,:,i);
    
    if i>1
        PHTelem(elementIndex).neighbor_left = indexMatrix(i-1);
    end
    if i<numElemU
        PHTelem(elementIndex).neighbor_right = indexMatrix(i+1);
    end    
end
[PHTelem,controlPts,dimBasis] = crossInsertIso1D(PHTelem,controlPts,1:numElements,dimBasis,p);
end