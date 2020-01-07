function [PHTelem,sizeBasis] = zipConforming3D(PHTelem,dimBasis,geometry)
% Connects two conforming patches by changing the nodesGlobal entry
% patchBoundaries format:
% patchA, patchB, faceA, faceB
% patchA should be 1
% face format: 1-front, 2-right, 3-back, 4-left, 5-down, 6-up

p = geometry.p;
q = geometry.q;
r = geometry.r;
patchBoundaries = geometry.patchBoundaries;

numBoundaries = size(patchBoundaries,1);
numPatches = length(PHTelem);

%create/set nodesGlobal entries in all patches to be equal to local nodes
%entries
for patchIndex = 1:numPatches
    for elemIndex = 1:length(PHTelem{patchIndex})
        PHTelem{patchIndex}(elemIndex).nodesGlobal = PHTelem{patchIndex}(elemIndex).nodes;
    end
end

curShift = dimBasis(1);
overlapCounter = 0;
patchesSeen = [];


for boundaryIndex = 1:numBoundaries
    %get the nodes on the boundary edge in patchA and patchB
    patchAList = patchBoundaries{boundaryIndex,1};
    patchB = patchBoundaries{boundaryIndex,2};
    faceAList = patchBoundaries{boundaryIndex,3};
    faceBList = patchBoundaries{boundaryIndex,4};
    nodesA = [];
    nodesB = [];
    patchesSeen = union(patchesSeen, patchB);
    for indexPatch=1:length(patchAList)
        
        patchA = patchAList(indexPatch);
        faceA = faceAList(indexPatch);
        faceB = faceBList(indexPatch);
        patchesSeen = union(patchesSeen, patchA);
        
        nodesPattern = zeros(1, dimBasis(patchB));
        nodesAtemp = sortFaceNodesElem(PHTelem{patchA},faceA,p,q,r);
        nodesBtemp = sortFaceNodesElem(PHTelem{patchB},faceB,p,q,r);
        
        nodesA = [nodesA, nodesAtemp];
        nodesB = [nodesB, nodesBtemp];
        if length(nodesA)~=length(nodesB)
            error('Non-conforming patches encountered. Aborting...')
        end
    end
    [nodesB,sI] = unique(nodesB);
    nodesA = nodesA(sI);
    
    curBdryNode = 0;
    for nodeIndex=1:length(nodesA)
        %shift the basis functions indices in nodesPattern
        prevBdryNode = curBdryNode;
        curBdryNode = nodesB(nodeIndex);
        nodesPattern(prevBdryNode+1:curBdryNode-1) = ((prevBdryNode+1):(curBdryNode-1)) + curShift;
        nodesPattern(curBdryNode) = nodesA(nodeIndex);
        curShift = curShift - 1;
    end
    %shift the indices after the last boundary node
    nodesPattern(curBdryNode+1:end) = ((curBdryNode+1):dimBasis(patchB)) + curShift;
    
    %update the nodesGlobal in patchB according to nodesPattern
    for elemIndex = 1:length(PHTelem{patchB})
        PHTelem{patchB}(elemIndex).nodesGlobal = nodesPattern(PHTelem{patchB}(elemIndex).nodes);
    end
    overlapCounter = overlapCounter + length(nodesA);
    curShift = sum(dimBasis(patchesSeen))-overlapCounter;
end
sizeBasis = sum(dimBasis) - overlapCounter;



