function [PHTelem,controlPts,Mesh,fieldData,markRef] = refineElemProjGradedIso1D(elemRef,PHTelem,controlPts,Mesh,fieldData,patchIndex)
% Refines the elements marked by elemRef also refine the neighbor patches if the level_neighbor<level_refined_quad
% Project the coefficients included in the field data
p = Mesh.p;
numElem = length(PHTelem); %contains the deactivated elements
numActiveElem = Mesh.numberElements; %only the active element
dimBasis = Mesh.dimBasis(patchIndex);
% This array has a marking of 0 for elements that
% are not refined and the ones that are not neighbouring to the elements that are refined.
markRef = zeros(length(PHTelem),1);  
% Check the refinement level for the neighbor quad and mark it for refinement if needed

% Sort elements by level to prevent coarse->fine level errors
% We must refine the coarsest marked elements first
levelList = zeros(1,numElem);
for i=1:numElem
    levelList(i) = PHTelem(i).level;
end

[~,sortedIndex] = sort(levelList);
elemList = 1:length(elemRef);
elemList = elemList(sortedIndex);
elemRef = elemRef(sortedIndex);
tempPHTelem = PHTelem;
for i=1:length(elemRef)
    if elemRef(i)
        curElem = elemList(i);
        if isempty(PHTelem(curElem).children)
            markRef(curElem) = 1;% The value of markRef for this element that is refined is marked 1
            east = tempPHTelem(curElem).neighbor_right;
            west = tempPHTelem(curElem).neighbor_left;
            neighbor = [east,west];
            % The neighbours of the element that has been refined is also 
            % marked since there is C1 continuity, so there is sharing of nodes
            % between neighbouring elements. MarkRef would be used as a reference
            % for the elements that have been changed in the calculation of 
            % the derivatives and also the transfer function
            markRef(neighbor) = 1; 
            [PHTelem,controlPts,dimBasis,fieldData] = crossInsertProjIso1D(PHTelem,controlPts,curElem,dimBasis,p,fieldData);
            % One new element created (1 deactivated parent and 2 children)
            numActiveElem = numActiveElem + 1;            
        else
            error('Something wrong')
        end
        
    end    
end
Mesh.dimBasis(patchIndex) = dimBasis;
Mesh.numberElements = numActiveElem;
