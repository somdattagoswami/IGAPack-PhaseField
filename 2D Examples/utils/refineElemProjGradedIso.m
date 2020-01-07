function [PHTelem,controlPts,dimBasis,fieldData,numberElements,markRef] = refineElemProjGradedIso(elemRef,PHTelem,controlPts,geometry,dimBasis,fieldData,numberElements)
% Refines the elements marked by elemRef also refine the neighbor patches if the level_neighbor<level_refined_quad
% Project the coefficients included in the field data

p = geometry.p;
q = geometry.q;

numElem = length(PHTelem);
keepChecking = 1;
markRef = zeros(length(PHTelem),1); % This array has a marking of 0 for elemnts that are not refined the ones that are not neighbouring to the elemnets that are refined.
% Check the refinement level for the neighbor quad and mark it for refinement if needed
extraElemCounter = 0;
while keepChecking
    keepChecking = 0;
    elemRefOrig = elemRef;
    for i=1:length(elemRefOrig)
        if elemRefOrig(i) && isempty(PHTelem(i).children)
            curLevel = PHTelem(i).level;
            if ~isempty(PHTelem(i).neighbor_down)
                neighborSouth = PHTelem(i).neighbor_down(1);
                if curLevel>PHTelem(neighborSouth).level
                    if elemRef(neighborSouth)==0
                        elemRef(neighborSouth) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
            if ~isempty(PHTelem(i).neighbor_left)
                neighborWest = PHTelem(i).neighbor_left(1);
                if curLevel>PHTelem(neighborWest).level
                    if elemRef(neighborWest)==0
                        elemRef(neighborWest) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
            if ~isempty(PHTelem(i).neighbor_up)
                neighborNorth = PHTelem(i).neighbor_up(1);
                if curLevel>PHTelem(neighborNorth).level
                    if elemRef(neighborNorth)==0
                        elemRef(neighborNorth) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
            
            if ~isempty(PHTelem(i).neighbor_right)
                neighborEast = PHTelem(i).neighbor_right(1);
                if curLevel>PHTelem(neighborEast).level
                    if elemRef(neighborEast)==0
                        elemRef(neighborEast) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
        end
    end
    disp(['Extra elements refined: ',num2str(extraElemCounter)]);
end


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
            north = tempPHTelem(curElem).neighbor_up;
            south = tempPHTelem(curElem).neighbor_down;
            east = tempPHTelem(curElem).neighbor_right;
            west = tempPHTelem(curElem).neighbor_left;
            if ~isempty(north)
                northeast = tempPHTelem(max(north)).neighbor_right;
                northwest = tempPHTelem(min(north)).neighbor_left;
            else
                northeast = [];
                northwest = [];
            end
            if ~isempty(south)
                southeast = tempPHTelem(max(south)).neighbor_right;
                southwest = tempPHTelem(min(south)).neighbor_left;
            else
                southeast = [];
                southwest = [];
            end
            neighbor = [east,west,north,south,northeast,northwest,southeast,southwest];
            markRef(neighbor) = 1; % The neighbours of the elemnent that has been refined is also marks since there is C1 continuity, so there is sharing of nodes between neighbouring elements. MarkRef would be used as a reference for the elemnts that have been changed incalculation of derivatives and also the transfer function
            [PHTelem, controlPts, dimBasis,fieldData] = crossInsertProjIso(PHTelem, controlPts,curElem,dimBasis,p,q,fieldData);
            % Three new elements created (1 deactivated parent and 4 children)
            numberElements = numberElements + 3;
            
        else
            warning('Something wrong')
            pause
        end
        
    end
end
