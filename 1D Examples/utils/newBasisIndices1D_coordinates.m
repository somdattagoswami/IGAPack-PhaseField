function [newBasisSet,dimBasis,newBasisIndex] = newBasisIndices1D_coordinates(elmIn,p,dimBasis)
%returns the entries of element nod after adding new basis vertices

alpha = floor((p-1)/2);
[west,center,east] = getSideIndices(p);
newBasisIndex = zeros(1,2);
%copy the entries in elmIn to each element in newBasisSet
newBasisSet = cell(1,2);
for i=1:2
    newBasisSet{i} = elmIn;
end

basisCounter = dimBasis;
%add the p-2*alpha-1 deleted basis functions
newBasisSet{1}(center) = elmIn(center);
%add the (alpha+1) basis functions around the center basis vertex
newBasisSet{1}(east) = basisCounter+1:basisCounter+1+alpha;
newBasisSet{2}(west) = basisCounter+1:basisCounter+1+alpha;
newBasisIndex(1,:) = basisCounter+1:basisCounter+1+alpha;
basisCounter = basisCounter+1+alpha;
%add p-2*alpha-1 basis functions in the center of the 2nd element
newBasisSet{2}(center) = basisCounter+1:basisCounter+p-2*alpha-1;
basisCounter = basisCounter+p-2*alpha-1;
%update dimension of basis space
dimBasis = basisCounter;

