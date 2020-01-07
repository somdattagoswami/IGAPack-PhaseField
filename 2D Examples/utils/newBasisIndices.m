function [newBasisSet,dimBasis] = newBasisIndices(newBasisVert,elmIn,p,q,dimBasis)
% Returns the entries of element nod after adding new basis vertices

[sw,south,se,west,center,east,nw,north,ne] = getCornerIndices(p,q);

alpha = floor((p-1)/2); beta = floor((q-1)/2);
% Copy the entries in elmIn to each element in newBasisSet
newBasisSet = cell(1,4);
for i=1:4
    newBasisSet{i} = elmIn;
end

basisCounter = dimBasis;
for j=1:length(newBasisVert)
    switch newBasisVert(j)
        case 1
            %add the 2*(p-3) deleted basis functions
            newBasisSet{1}(south) = elmIn(south);
            %add 4 more basis function near the vertex
            newBasisSet{1}(se) = basisCounter+1:basisCounter+(alpha+1)^2;
            newBasisSet{2}(sw) = basisCounter+1:basisCounter+(alpha+1)^2;
            %add 2*(p-3) more basis function in the middle of the south edge
            newBasisSet{2}(south) = basisCounter+(beta+1)^2+1:basisCounter+(alpha+1)^2+(beta+1)*(p-2*alpha-1);
            %increment basis counter
            basisCounter = basisCounter+(alpha+1)^2+(alpha+1)*(p-2*alpha-1);
        case 2
            %add the 2*(q-3) deleted basis functions
            newBasisSet{2}(east) = elmIn(east);
            %add 4 more basis function near the vertex
            newBasisSet{2}(ne) = basisCounter+1:basisCounter+(beta+1)^2;
            newBasisSet{4}(se) = basisCounter+1:basisCounter+(beta+1)^2;
            %add 2*(q-3) more basis function in the middle of the east edge
            newBasisSet{4}(east) = basisCounter+(beta+1)^2+1:basisCounter+(beta+1)^2+(alpha+1)*(q-2*beta-1);
            %increment basis counter
            basisCounter = basisCounter+(beta+1)^2+(alpha+1)*(q-2*beta-1);
        case 3
            %add the 2*(p-3) deleted basis functions
            newBasisSet{3}(north) = elmIn(north);
            %add 4 more basis function near the vertex
            newBasisSet{3}(ne) = basisCounter+1:basisCounter+(alpha+1)^2;
            newBasisSet{4}(nw) = basisCounter+1:basisCounter+(alpha+1)^2;
            %add 2*(p-3) more basis function in the middle of the north edge
            newBasisSet{4}(north) = basisCounter+(alpha+1)^2+1:basisCounter+(alpha+1)^2+(beta+1)*(p-2*alpha-1);
            %increment basis counter
            basisCounter = basisCounter+(alpha+1)^2+(beta+1)*(p-2*alpha-1);
        case 4
            %add the 2*(q-3) deleted basis functions
            newBasisSet{1}(west) = elmIn(west);
            %add 4 more basis function near the vertex
            newBasisSet{1}(nw) = basisCounter+1:basisCounter+(beta+1)^2;
            newBasisSet{3}(sw) = basisCounter+1:basisCounter+(beta+1)^2;
            %add 2*(q-3) more basis function in the middle of the west edge
            newBasisSet{3}(west) = basisCounter+(beta+1)^2+1:basisCounter+(beta+1)^2+(alpha+1)*(q-2*beta-1);
            %increment basis counter
            basisCounter = basisCounter+(alpha+1)^2+(alpha+1)*(q-2*beta-1);
        case 5
            %add (p-3)*(q-3) deleted basis functions
            newBasisSet{1}(center) = elmIn(center);
            %add the other new basis functions
            newBasisSet{1}(east) = basisCounter+1:basisCounter+(beta+1)*(q-2*beta-1);
            newBasisSet{2}(west) = basisCounter+1:basisCounter+(beta+1)*(q-2*beta-1);
            basisCounter = basisCounter + (beta+1)*(q-2*beta-1);
            newBasisSet{2}(center) = basisCounter+1:basisCounter+(p-2*alpha-1)*(q-2*beta-1);
            basisCounter = basisCounter + (p-2*alpha-1)*(q-2*beta-1);
            newBasisSet{1}(north) = basisCounter+1:basisCounter+(alpha+1)*(p-2*alpha-1);
            newBasisSet{3}(south) = basisCounter+1:basisCounter+(alpha+1)*(p-2*alpha-1);
            basisCounter = basisCounter + (alpha+1)*(p-2*alpha-1);
            newBasisSet{1}(ne) = basisCounter+1:basisCounter+(alpha+1)^2;
            newBasisSet{2}(nw) = basisCounter+1:basisCounter+(alpha+1)^2;
            newBasisSet{3}(se) = basisCounter+1:basisCounter+(alpha+1)^2;
            newBasisSet{4}(sw) = basisCounter+1:basisCounter+(alpha+1)^2;
            basisCounter = basisCounter + (alpha+1)^2;
            newBasisSet{3}(center) = basisCounter+1:basisCounter+(p-2*alpha-1)*(q-2*beta-1);
            basisCounter = basisCounter + (p-2*alpha-1)*(q-2*beta-1);
            newBasisSet{3}(east) = basisCounter+1:basisCounter+(beta+1)*(q-2*beta-1);
            newBasisSet{4}(west) = basisCounter+1:basisCounter+(beta+1)*(q-2*beta-1);
            basisCounter = basisCounter + (beta+1)*(q-2*beta-1);
            newBasisSet{4}(center) = basisCounter+1:basisCounter+(p-2*alpha-1)*(q-2*beta-1);
            basisCounter = basisCounter + (p-2*alpha-1)*(q-2*beta-1);
            newBasisSet{2}(north) = basisCounter+1:basisCounter+ (alpha+1)*(p-2*alpha-1);
            newBasisSet{4}(south) = basisCounter+1:basisCounter+ (alpha+1)*(p-2*alpha-1);
            basisCounter = basisCounter + (alpha+1)*(p-2*alpha-1);
            
    end
end

% Update dimension of basis space
dimBasis = basisCounter;
end

