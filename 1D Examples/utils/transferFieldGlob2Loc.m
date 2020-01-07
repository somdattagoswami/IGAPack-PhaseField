function [fieldDataPatch] = transferFieldGlob2Loc(PHTelem,dimBasis,fieldData)
%transfer the fieldData inidices from the global array to the local patch
%structure

numPatches = length(PHTelem);
fieldDataPatch = cell(1, numPatches);

for iPatch = 1:numPatches
    fieldDataPatch{iPatch} = zeros(dimBasis(iPatch),size(fieldData,2));
    for iElem = 1:length(PHTelem{iPatch})
        if isempty(PHTelem{iPatch}(iElem).children)
            globIndex = PHTelem{iPatch}(iElem).nodesGlobal;
            localIndex = PHTelem{iPatch}(iElem).nodes;
            fieldDataPatch{iPatch}(localIndex,:) = fieldData(globIndex,:);
        end
    end
end