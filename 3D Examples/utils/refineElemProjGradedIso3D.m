function [PHTelem,controlPts,dimBasis,fieldData,numberElements] = refineElemProjGradedIso3D(elemRef,PHTelem,controlPts,geometry,dimBasis,fieldData,numberElements)
% Refines the elements marked by elemRef also refine the neighbor patches if the level_neighbor<level_refined_quad
% Project the coefficients included in the field data
p = geometry.p;
q = geometry.q;
r = geometry.r;
numElem = length(PHTelem);
keepChecking = 1;
% Check the refinement level for the neighbor quad and mark it for refinement if needed
extraElemCounter = 0;
while keepChecking
    keepChecking = 0;
    elemRefOrig = elemRef;
    for i=1:length(elemRefOrig)
        if elemRefOrig(i) && isempty(PHTelem(i).children)
            curLevel = PHTelem(i).level;
            if ~isempty(PHTelem(i).neighbor_left)
                neighborLeft = PHTelem(i).neighbor_left(1);
                if curLevel>PHTelem(neighborLeft).level
                    if elemRef(neighborLeft)==0
                        elemRef(neighborLeft) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
            if ~isempty(PHTelem(i).neighbor_right)
                neighborRight = PHTelem(i).neighbor_right(1);
                if curLevel>PHTelem(neighborRight).level
                    if elemRef(neighborRight)==0
                        elemRef(neighborRight) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
            if ~isempty(PHTelem(i).neighbor_down)
                neighborDown = PHTelem(i).neighbor_down(1);
                if curLevel > PHTelem(neighborDown).level
                    if elemRef(neighborDown)==0
                        elemRef(neighborDown) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end  
            if ~isempty(PHTelem(i).neighbor_up)
                neighborUp = PHTelem(i).neighbor_up(1);
                if curLevel>PHTelem(neighborUp).level
                    if elemRef(neighborUp)==0
                        elemRef(neighborUp) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
            if ~isempty(PHTelem(i).neighbor_front)
                neighborFront = PHTelem(i).neighbor_front(1);
                if curLevel > PHTelem(neighborFront).level
                    if elemRef(neighborFront)==0
                        elemRef(neighborFront) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
            if ~isempty(PHTelem(i).neighbor_back)
                neighborBack = PHTelem(i).neighbor_back(1);
                if curLevel>PHTelem(neighborBack).level
                    if elemRef(neighborBack)==0
                        elemRef(neighborBack) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
            if ~isempty(PHTelem(i).neighbor_up_left)
                neighborUpLeft = PHTelem(i).neighbor_up_left(1);
                if curLevel>PHTelem(neighborUpLeft).level
                    if elemRef(neighborUpLeft)==0
                        elemRef(neighborUpLeft) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
            if ~isempty(PHTelem(i).neighbor_down_left)
                neighborDownLeft = PHTelem(i).neighbor_down_left(1);
                if curLevel>PHTelem(neighborDownLeft).level
                    if elemRef(neighborDownLeft)==0
                        elemRef(neighborDownLeft) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
            if ~isempty(PHTelem(i).neighbor_up_right)
                neighborUpRight = PHTelem(i).neighbor_up_right(1);
                if curLevel>PHTelem(neighborUpRight).level
                    if elemRef(neighborUpRight)==0
                        elemRef(neighborUpRight) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
            if ~isempty(PHTelem(i).neighbor_down_right)
                neighborDownRight = PHTelem(i).neighbor_down_right(1);
                if curLevel>PHTelem(neighborDownRight).level
                    if elemRef(neighborDownRight)==0
                        elemRef(neighborDownRight) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
            if ~isempty(PHTelem(i).neighbor_up_front)
                neighborUpFront = PHTelem(i).neighbor_up_front(1);
                if curLevel>PHTelem(neighborUpFront).level
                    if elemRef(neighborUpFront)==0
                        elemRef(neighborUpFront) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
            if ~isempty(PHTelem(i).neighbor_down_front)
                neighborDownFront = PHTelem(i).neighbor_down_front(1);
                if curLevel>PHTelem(neighborDownFront).level
                    if elemRef(neighborDownFront)==0
                        elemRef(neighborDownFront) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
            if ~isempty(PHTelem(i).neighbor_up_back)
                neighborUpBack = PHTelem(i).neighbor_up_back(1);
                if curLevel>PHTelem(neighborUpBack).level
                    if elemRef(neighborUpBack)==0
                        elemRef(neighborUpBack) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
            if ~isempty(PHTelem(i).neighbor_down_back)
                neighborDownBack = PHTelem(i).neighbor_down_back(1);
                if curLevel>PHTelem(neighborDownBack).level
                    if elemRef(neighborDownBack)==0
                        elemRef(neighborDownBack) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end    
            if ~isempty(PHTelem(i).neighbor_left_front)
                neighborLeftFront = PHTelem(i).neighbor_left_front(1);
                if curLevel>PHTelem(neighborLeftFront).level
                    if elemRef(neighborLeftFront)==0
                        elemRef(neighborLeftFront) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
            if ~isempty(PHTelem(i).neighbor_right_front)
                neighborRightFront = PHTelem(i).neighbor_right_front(1);
                if curLevel>PHTelem(neighborRightFront).level
                    if elemRef(neighborRightFront)==0
                        elemRef(neighborRightFront) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
            if ~isempty(PHTelem(i).neighbor_left_back)
                neighborLeftBack = PHTelem(i).neighbor_left_back(1);
                if curLevel>PHTelem(neighborLeftBack).level
                    if elemRef(neighborLeftBack)==0
                        elemRef(neighborLeftBack) = 1;
                        extraElemCounter = extraElemCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
            if ~isempty(PHTelem(i).neighbor_right_back)
                neighborRightBack = PHTelem(i).neighbor_right_back(1);
                if curLevel>PHTelem(neighborRightBack).level
                    if elemRef(neighborRightBack)==0
                        elemRef(neighborRightBack) = 1;
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
for i=1:length(elemRef)
    if elemRef(i)
        curElem = elemList(i);
        if isempty(PHTelem(curElem).children)
            [PHTelem,controlPts,dimBasis,fieldData] = crossInsertProjIso3D(PHTelem,controlPts,curElem,dimBasis,p,q,r,fieldData);
            % Seven new elements created (1 deactivated parent and 8 children)
            numberElements = numberElements + 7;            
        else
            warning('Something wrong')
            pause
        end
        
    end
end
