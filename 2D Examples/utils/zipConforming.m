function [PHTelem,sizeBasis] = zipConforming(PHTelem,dimBasis,geometry)
%connects two conforming patches by changing the nodesGlobal entry
%patchBoundaries format:
% patchA, patchB, edgeA, edgeB, flagRepeat
% patchA should be 1
% flagRepeat = 0 --> shift node indices, flagRepeat = 1 -->
%edge format: 1-down, 2-right, 3-up, 4-left
p = geometry.p;
q = geometry.q;
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
    edgeAList = patchBoundaries{boundaryIndex,3};
    edgeBList = patchBoundaries{boundaryIndex,4};
    nodesPattern = zeros(1, dimBasis(patchB));
    nodesA = [];
    nodesB = [];
    patchesSeen = union(patchesSeen, patchB);
    
    for indexPatch=1:length(patchAList)
        patchA = patchAList(indexPatch);
        edgeA = edgeAList(indexPatch);
        edgeB = edgeBList(indexPatch);
        patchesSeen = union(patchesSeen, patchA);
        
        nodesAcell = sortEdgeNodesElem( PHTelem{patchA}, edgeA, p, q );
        nodesBcell = sortEdgeNodesElem( PHTelem{patchB}, edgeB, p, q );
        
        nodesA = [nodesA, nodesAcell{1}];
        nodesB = [nodesB, nodesBcell{1}];
        
        if length(nodesA)~=length(nodesB)
            error('Non-conforming patches encountered. Aborting...')
        end
    end
    [nodesB,sI] = sort(nodesB);
    nodesA = nodesA(sI);

    curBdryNode = 0;
    for nodeIndex=1:length(nodesA)
        %shift the basis functions indices in nodesPattern
        prevBdryNode = curBdryNode;
        curBdryNode = nodesB(nodeIndex);
        nodesPattern(prevBdryNode+1:curBdryNode-1) = ((prevBdryNode+1):(curBdryNode-1)) + curShift;
        nodesPattern(curBdryNode) = nodesA(nodeIndex);
        if prevBdryNode<curBdryNode
            curShift = curShift - 1;
        end
    end
    %shift the indices after the last boundary node
    nodesPattern(curBdryNode+1:end) = ((curBdryNode+1):dimBasis(patchB)) + curShift;
    %update the nodesGlobal in patchB according to nodesPattern
    for elemIndex = 1:length(PHTelem{patchB})
        PHTelem{patchB}(elemIndex).nodesGlobal = nodesPattern(PHTelem{patchB}(elemIndex).nodes);
    end    
    overlapCounter = overlapCounter + length(unique(nodesA));
    
    curShift = sum(dimBasis(patchesSeen))-overlapCounter;    
    
end
sizeBasis = sum(dimBasis) - overlapCounter;



