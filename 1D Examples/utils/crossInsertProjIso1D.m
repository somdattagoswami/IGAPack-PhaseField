function [PHTelem,controlPts,dimBasis,fieldData] = crossInsertProjIso1D(PHTelem,controlPts,elm_index,dimBasis,p,fieldData)
%inserts cross in PHTelem at elm_index
%dimBasis: dimension of the PHT space
%projects field data on the new mesh

coordinates = [controlPts,fieldData];

for e=elm_index
    lastElem = length(PHTelem);    
    
    % Update parent element
    PHTelem(e).children = lastElem+1:lastElem+2;
    umin = PHTelem(e).vertex(1);
    umax = PHTelem(e).vertex(2);
    ucenter = (umin+umax)/2;
    numEnt = size(PHTelem(e).C,1);
    nodes = PHTelem(e).nodes(1:numEnt);
    cpts = coordinates(nodes(1:numEnt),:);
    
    geomInfo = getGeometricInformation1D(PHTelem(e),cpts,p,umin,umax);
    %geomInfo{1:end}
    
    DeltaU1 = (ucenter - umin)/3;
    DeltaU2 = (umax - ucenter)/3;
    alpha=1/(DeltaU1+DeltaU2);
    lambda=alpha*DeltaU1;    
    b =[(1-lambda),lambda; -alpha,alpha];
    rhs = geomInfo;
    pcpts = b\rhs;
    
    % Add children elements and their child neighbors
    % Left child
    PHTelem(lastElem+1).parent = e;
    PHTelem(lastElem+1).children = [];
    PHTelem(lastElem+1).vertex = [umin, (umin+umax)/2];
    PHTelem(lastElem+1).level = PHTelem(e).level+1;   
    PHTelem(lastElem+1).neighbor_right = lastElem+2;    
    
    %right child
    PHTelem(lastElem+2).parent = e;
    PHTelem(lastElem+2).children = [];
    PHTelem(lastElem+2).vertex = [(umin+umax)/2, umax];
    PHTelem(lastElem+2).level = PHTelem(e).level+1;
    PHTelem(lastElem+2).neighbor_left = lastElem+1;
   
    %add the neighbors outside the refined element
    PHTelem(lastElem+2).neighbor_right = PHTelem(e).neighbor_right;
    PHTelem(lastElem+1).neighbor_left = PHTelem(e).neighbor_left;
    
    %calculate the new Bezier extraction operators and element_nod indices
    %of the children elements
    C_temp = PHTelem(e).C;
    elmIn = PHTelem(e).nodes;
    
    leftNeighbor = PHTelem(e).neighbor_left;
    if ~isempty(leftNeighbor)
        knotUl = PHTelem(leftNeighbor).vertex(1);
    else
        knotUl = 0;
    end
    rightNeighbor = PHTelem(e).neighbor_right;
    if ~isempty(rightNeighbor)        
        knotUr = PHTelem(rightNeighbor).vertex(2);
    else
        knotUr = 1;
    end
    [Ce1,Ce2,in1,in2,dimBasis,newBasisIndex] = deCasteljau1dai_coordinates(C_temp,umin,umax,p,knotUl,knotUr,elmIn,dimBasis);    
    coordinates(newBasisIndex(1,:),:) = zeros(2,size(coordinates,2));
    coordinates(newBasisIndex(1,:),1) = pcpts(:,1)./pcpts(:,2);
    coordinates(newBasisIndex(1,:),2) = pcpts(:,2);
    coordinates(newBasisIndex(1,:),3:end) = pcpts(:,3:end);
    
     PHTelem(lastElem+1).C = Ce1;
    PHTelem(lastElem+2).C = Ce2;
    
    PHTelem(lastElem+1).nodes = in1;
    PHTelem(lastElem+2).nodes = in2;
                   
    %update the neighbors of the neighbors with self
    for ichild = 1:2        
        if (length(PHTelem(lastElem+ichild).neighbor_right)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_right).neighbor_left = setdiff(unique([PHTelem(PHTelem(lastElem+ichild).neighbor_right).neighbor_left, lastElem+ichild]),e);
        end      
        
        if (length(PHTelem(lastElem+ichild).neighbor_left)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_left).neighbor_right = setdiff(unique([PHTelem(PHTelem(lastElem+ichild).neighbor_left).neighbor_right, lastElem+ichild]),e);
        end
    end                    
end
%TODO: Change if using rational functions
controlPts = coordinates(:,1:2);
fieldData = coordinates(:,3:end);