function[PHTelem,sizeBasis] = zipConforming1D(PHTelem,dimBasis)

numPatches = length(PHTelem);
for patchIndex = 1:numPatches
    for elemIndex = 1:length(PHTelem{patchIndex})
        PHTelem{patchIndex}(elemIndex).nodesGlobal = PHTelem{patchIndex}(elemIndex).nodes;
    end
end

overlapCounter = 0;
sizeBasis = sum(dimBasis) - overlapCounter;
end