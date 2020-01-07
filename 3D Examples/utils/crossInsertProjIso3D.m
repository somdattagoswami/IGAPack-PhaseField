function [PHTelem,controlPts,dimBasis,fieldData] = crossInsertProjIso3D(PHTelem,controlPts,elm_index,dimBasis,p,q,r,fieldData)
%inserts cross in PHTelem at elm_index
%dimBasis: dimension of the PHT space
%projects field data on the new mesh

coordinates = [controlPts,fieldData];

% Define corner indices
reg = getRegionIndices3D(p,q,r);

for e = elm_index
    lastElem = length(PHTelem);
    % Collect information about new basis vertices and T junctions
    [newBasisVert,~,RTjunct,knotUl,knotUr,knotVd,knotVu,knotWb,knotWt] = checkNeighbors3D(PHTelem,e);
    
    % Update parent element
    PHTelem(e).children = lastElem+1:lastElem+8;
    
    xmin = PHTelem(e).vertex(1);
    xmax = PHTelem(e).vertex(4);
    ymin = PHTelem(e).vertex(2);
    ymax = PHTelem(e).vertex(5);
    zmin = PHTelem(e).vertex(3);
    zmax = PHTelem(e).vertex(6);
    xcenter = (xmin+xmax)/2;
    ycenter = (ymin+ymax)/2;
    zcenter = (zmin+zmax)/2;
    nument = size(PHTelem(e).C,1);
    nodes = PHTelem(e).nodes(1:nument);
    cpts = coordinates(nodes,:);
    
    geomInfo = getGeometricInformation3D(PHTelem(e),cpts,[newBasisVert,RTjunct,7],p,q,r,xmin,xmax,ymin,ymax,zmin,zmax);
        
    indexCounter = 0;
    cpts = cell(length([newBasisVert,RTjunct,7]),1);
    
    for indexNewVert = [newBasisVert,RTjunct,7]
        indexCounter = indexCounter + 1;
        
        switch indexNewVert
            
            case 1
                
                deltaU1=xcenter-xmin;
                deltaU2=xmax-xcenter;
                deltaV1=ymin-knotVd;
                deltaV2=ycenter-ymin;
                deltaW1=zcenter-zmin;
                deltaW2=zmax-zcenter;
                
            case 2
                
                deltaU1=xmax-xcenter;
                deltaU2=knotUr-xmax;
                deltaV1=ycenter-ymin;
                deltaV2=ymax-ycenter;
                deltaW1=zcenter-zmin;
                deltaW2=zmax-zcenter;
                
            case 3
                
                deltaU1=xcenter-xmin;
                deltaU2=xmax-xcenter;
                deltaV1=ymax-ycenter;
                deltaV2=knotVu-ymax;
                deltaW1=zcenter-zmin;
                deltaW2=zmax-zcenter;
                
            case 4
                
                deltaU1=xmin-knotUl;
                deltaU2=xcenter-xmin;
                deltaV1=ycenter-ymin;
                deltaV2=ymax-ycenter;
                deltaW1=zcenter-zmin;
                deltaW2=zmax-zcenter;
                
            case 5
                
                deltaU1=xcenter-xmin;
                deltaU2=xmax-xcenter;
                deltaV1=ycenter-ymin;
                deltaV2=ymax-ycenter;
                deltaW1=zmin-knotWb;
                deltaW2=zcenter-zmin;
                
            case 6
                
                deltaU1=xcenter-xmin;
                deltaU2=xmax-xcenter;
                deltaV1=ycenter-ymin;
                deltaV2=ymax-ycenter;
                deltaW1=zmax-zcenter;
                deltaW2=knotWt-zmax;
                
            case 7
                
                deltaU1=xcenter-xmin;
                deltaU2=xmax-xcenter;
                deltaV1=ycenter-ymin;
                deltaV2=ymax-ycenter;
                deltaW1=zcenter-zmin;
                deltaW2=zmax-zcenter;
                
            case 8
                
                deltaU1=xmin-knotUl;
                deltaU2=xcenter-xmin;
                deltaV1=ycenter-ymin;
                deltaV2=ymax-ycenter;
                deltaW1=zmax-zcenter;
                deltaW2=knotWt-zmax;
                
            case 9
                
                deltaU1=xmin-knotUl;
                deltaU2=xcenter-xmin;
                deltaV1=ycenter-ymin;
                deltaV2=ymax-ycenter;
                deltaW1=zmin-knotWb;
                deltaW2=zcenter-zmin;
                
            case 10
                
                deltaU1=xmax-xcenter;
                deltaU2=knotUr-xmax;
                deltaV1=ycenter-ymin;
                deltaV2=ymax-ycenter;
                deltaW1=zmax-zcenter;
                deltaW2=knotWt-zmax;
                
            case 11
                
                deltaU1=xmax-xcenter;
                deltaU2=knotUr-xmax;
                deltaV1=ycenter-ymin;
                deltaV2=ymax-ycenter;
                deltaW1=zmin-knotWb;
                deltaW2=zcenter-zmin;
                
            case 12
                
                deltaU1=xcenter-xmin;
                deltaU2=xmax-xcenter;
                deltaV1=ymin-knotVd;
                deltaV2=ycenter-ymin;
                deltaW1=zmax-zcenter;
                deltaW2=knotWt-zmax;
                
            case 13
                
                
                deltaU1=xcenter-xmin;
                deltaU2=xmax-xcenter;
                deltaV1=ymin-knotVd;
                deltaV2=ycenter-ymin;
                deltaW1=zmin-knotWb;
                deltaW2=zcenter-zmin;
                
            case 14
                
                deltaU1=xcenter-xmin;
                deltaU2=xmax-xcenter;
                deltaV1=ymax-ycenter;
                deltaV2=knotVu-ymax;
                deltaW1=zmax-zcenter;
                deltaW2=knotWt-zmax;
                
            case 15
                
                deltaU1=xcenter-xmin;
                deltaU2=xmax-xcenter;
                deltaV1=ymax-ycenter;
                deltaV2=knotVu-ymax;
                deltaW1=zmin-knotWb;
                deltaW2=zcenter-zmin;
                
            case 16
                
                deltaU1=xmin-knotUl;
                deltaU2=xcenter-xmin;
                deltaV1=ymin-knotVd;
                deltaV2=ycenter-ymin;
                deltaW1=zcenter-zmin;
                deltaW2=zmax-zcenter;
                
            case 17
                
                deltaU1=xmax-xcenter;
                deltaU2=knotUr-xmax;
                deltaV1=ymin-knotVd;
                deltaV2=ycenter-ymin;
                deltaW1=zcenter-zmin;
                deltaW2=zmax-zcenter;
                
            case 18
                
                deltaU1=xmin-knotUl;
                deltaU2=xcenter-xmin;
                deltaV1=ymax-ycenter;
                deltaV2=knotVu-ymax;
                deltaW1=zcenter-zmin;
                deltaW2=zmax-zcenter;
                
            case 19
                
                deltaU1=xmax-xcenter;
                deltaU2=knotUr-xmax;
                deltaV1=ymax-ycenter;
                deltaV2=knotVu-ymax;
                deltaW1=zcenter-zmin;
                deltaW2=zmax-zcenter;
        end
        
        deltaV1=deltaV1/3;
        deltaV2=deltaV2/3;
        deltaU1=deltaU1/3;
        deltaU2=deltaU2/3;
        deltaW1=deltaW1/3;
        deltaW2=deltaW2/3;
        
        alpha=1/(deltaU1+deltaU2);
        beta=1/(deltaV1+deltaV2);
        gamma = 1/(deltaW1+deltaW2);
        lambda=alpha*deltaU1;
        mu=beta*deltaV1;
        nu=gamma*deltaW1;
        
        blockN = [(1-lambda),-alpha;lambda,alpha]';
        blockM = [(1-mu),-beta; mu,beta]';
        blockP = [(1-nu),-gamma; nu,gamma]';
        
        b = kron(blockP,kron(blockM,blockN));
        b(4:5,:) = b([5,4],:);
                
        rhs = geomInfo{indexCounter};
        cpts{indexCounter} = b\rhs;
        
        
    end
    
    %add children elements and their easy neighbors
    %lower-front-left child
    PHTelem(lastElem+1).parent = e;
    PHTelem(lastElem+1).children = [];
    PHTelem(lastElem+1).vertex = [xmin, ymin, zmin, (xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2];
    PHTelem(lastElem+1).level = PHTelem(e).level+1;
    PHTelem(lastElem+1).neighbor_up = lastElem+5;
    PHTelem(lastElem+1).neighbor_right = lastElem+2;
    PHTelem(lastElem+1).neighbor_back = lastElem+3;
    PHTelem(lastElem+1).neighbor_up_right = lastElem+6;
    PHTelem(lastElem+1).neighbor_right_back = lastElem+4;
    PHTelem(lastElem+1).neighbor_up_back  = lastElem+7;
    PHTelem(lastElem+1).level = PHTelem(e).level+1;
    
    %lower-front-right child
    PHTelem(lastElem+2).parent = e;
    PHTelem(lastElem+2).children = [];
    PHTelem(lastElem+2).vertex = [(xmin+xmax)/2, ymin, zmin, xmax, (ymin+ymax)/2, (zmin+zmax)/2];
    PHTelem(lastElem+2).level = PHTelem(e).level+1;
    PHTelem(lastElem+2).neighbor_up = lastElem+6;
    PHTelem(lastElem+2).neighbor_left = lastElem+1;
    PHTelem(lastElem+2).neighbor_back = lastElem+4;
    PHTelem(lastElem+2).neighbor_up_left = lastElem+5;
    PHTelem(lastElem+2).neighbor_left_back = lastElem+3;
    PHTelem(lastElem+2).neighbor_up_back  = lastElem+8;
    PHTelem(lastElem+2).level = PHTelem(e).level+1;
    
    %lower-back-left child
    PHTelem(lastElem+3).parent = e;
    PHTelem(lastElem+3).children = [];
    PHTelem(lastElem+3).vertex = [xmin, (ymin+ymax)/2, zmin, (xmin+xmax)/2, ymax, (zmin+zmax)/2];
    PHTelem(lastElem+3).level = PHTelem(e).level+1;
    PHTelem(lastElem+3).neighbor_up = lastElem+7;
    PHTelem(lastElem+3).neighbor_right = lastElem+4;
    PHTelem(lastElem+3).neighbor_front = lastElem+1;
    PHTelem(lastElem+3).neighbor_up_front = lastElem+5;
    PHTelem(lastElem+3).neighbor_right_front = lastElem+2;
    PHTelem(lastElem+3).neighbor_up_right  = lastElem+8;
    PHTelem(lastElem+3).level = PHTelem(e).level+1;
        
    %lower-back-right child
    PHTelem(lastElem+4).parent = e;
    PHTelem(lastElem+4).children = [];
    PHTelem(lastElem+4).vertex = [(xmin+xmax)/2, (ymin+ymax)/2, zmin, xmax, ymax, (zmin+zmax)/2];
    PHTelem(lastElem+4).level = PHTelem(e).level+1;
    PHTelem(lastElem+4).neighbor_up = lastElem+8;
    PHTelem(lastElem+4).neighbor_left = lastElem+3;
    PHTelem(lastElem+4).neighbor_front = lastElem+2;
    PHTelem(lastElem+4).neighbor_up_front = lastElem+6;
    PHTelem(lastElem+4).neighbor_left_front = lastElem+1;
    PHTelem(lastElem+4).neighbor_up_left  = lastElem+7;
    PHTelem(lastElem+4).level = PHTelem(e).level+1;
        
    %upper-front-left child
    PHTelem(lastElem+5).parent = e;
    PHTelem(lastElem+5).children = [];
    PHTelem(lastElem+5).vertex = [xmin, ymin, (zmin+zmax)/2, (xmin+xmax)/2, (ymin+ymax)/2, zmax];
    PHTelem(lastElem+5).level = PHTelem(e).level+1;
    PHTelem(lastElem+5).neighbor_down = lastElem+1;
    PHTelem(lastElem+5).neighbor_right = lastElem+6;
    PHTelem(lastElem+5).neighbor_back = lastElem+7;
    PHTelem(lastElem+5).neighbor_down_right = lastElem+2;
    PHTelem(lastElem+5).neighbor_right_back = lastElem+8;
    PHTelem(lastElem+5).neighbor_down_back  = lastElem+3;
    PHTelem(lastElem+5).level = PHTelem(e).level+1;    
    PHTelem(lastElem+5).adjVertexlist=[];
    
    %upper-front-right child
    PHTelem(lastElem+6).parent = e;
    PHTelem(lastElem+6).children = [];
    PHTelem(lastElem+6).vertex = [(xmin+xmax)/2, ymin, (zmin+zmax)/2, xmax, (ymin+ymax)/2, zmax];
    PHTelem(lastElem+6).level = PHTelem(e).level+1;
    PHTelem(lastElem+6).neighbor_down = lastElem+2;
    PHTelem(lastElem+6).neighbor_left = lastElem+5;
    PHTelem(lastElem+6).neighbor_back = lastElem+8;
    PHTelem(lastElem+6).neighbor_down_left = lastElem+1;
    PHTelem(lastElem+6).neighbor_left_back = lastElem+7;
    PHTelem(lastElem+6).neighbor_down_back  = lastElem+4;
    PHTelem(lastElem+6).level = PHTelem(e).level+1;
    
    %upper-back-left child
    PHTelem(lastElem+7).parent = e;
    PHTelem(lastElem+7).children = [];
    PHTelem(lastElem+7).vertex = [xmin, (ymin+ymax)/2, (zmin+zmax)/2, (xmin+xmax)/2, ymax, zmax];
    PHTelem(lastElem+7).level = PHTelem(e).level+1;
    PHTelem(lastElem+7).neighbor_down = lastElem+3;
    PHTelem(lastElem+7).neighbor_right = lastElem+8;
    PHTelem(lastElem+7).neighbor_front = lastElem+5;
    PHTelem(lastElem+7).neighbor_down_front = lastElem+1;
    PHTelem(lastElem+7).neighbor_right_front = lastElem+6;
    PHTelem(lastElem+7).neighbor_down_right  = lastElem+4;
    PHTelem(lastElem+7).level = PHTelem(e).level+1;
        
    %upper-back-right child
    PHTelem(lastElem+8).parent = e;
    PHTelem(lastElem+8).children = [];
    PHTelem(lastElem+8).vertex = [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2, xmax, ymax, zmax];
    PHTelem(lastElem+8).level = PHTelem(e).level+1;
    PHTelem(lastElem+8).neighbor_down = lastElem+4;
    PHTelem(lastElem+8).neighbor_left = lastElem+7;
    PHTelem(lastElem+8).neighbor_front = lastElem+6;
    PHTelem(lastElem+8).neighbor_down_front = lastElem+2;
    PHTelem(lastElem+8).neighbor_left_front = lastElem+5;
    PHTelem(lastElem+8).neighbor_down_left  = lastElem+3;
    PHTelem(lastElem+8).level = PHTelem(e).level+1; 
    
   
   %add the neighbors outside the refined element, taking into account
    %T-junctions
    
    if (length(PHTelem(e).neighbor_down)==4)
        PHTelem(lastElem+1).neighbor_down = PHTelem(e).neighbor_down(1);
        PHTelem(lastElem+2).neighbor_down = PHTelem(e).neighbor_down(2);
        PHTelem(lastElem+3).neighbor_down = PHTelem(e).neighbor_down(3);
        PHTelem(lastElem+4).neighbor_down = PHTelem(e).neighbor_down(4);
        PHTelem(lastElem+1).neighbor_down_right = PHTelem(e).neighbor_down(2);
        PHTelem(lastElem+1).neighbor_down_back = PHTelem(e).neighbor_down(3);
        PHTelem(lastElem+2).neighbor_down_left = PHTelem(e).neighbor_down(1);
        PHTelem(lastElem+2).neighbor_down_back = PHTelem(e).neighbor_down(4);
        PHTelem(lastElem+3).neighbor_down_right = PHTelem(e).neighbor_down(4);
        PHTelem(lastElem+3).neighbor_down_front = PHTelem(e).neighbor_down(1);
        PHTelem(lastElem+4).neighbor_down_left = PHTelem(e).neighbor_down(3);
        PHTelem(lastElem+4).neighbor_down_front = PHTelem(e).neighbor_down(2);
    else
        PHTelem(lastElem+1).neighbor_down = PHTelem(e).neighbor_down;
        PHTelem(lastElem+2).neighbor_down = PHTelem(e).neighbor_down;
        PHTelem(lastElem+3).neighbor_down = PHTelem(e).neighbor_down;
        PHTelem(lastElem+4).neighbor_down = PHTelem(e).neighbor_down;
    end
    
    if (length(PHTelem(e).neighbor_right)==4)
        PHTelem(lastElem+2).neighbor_right = PHTelem(e).neighbor_right(1);
        PHTelem(lastElem+4).neighbor_right = PHTelem(e).neighbor_right(2);
        PHTelem(lastElem+6).neighbor_right = PHTelem(e).neighbor_right(3);
        PHTelem(lastElem+8).neighbor_right = PHTelem(e).neighbor_right(4);
        PHTelem(lastElem+2).neighbor_up_right = PHTelem(e).neighbor_right(3);
        PHTelem(lastElem+4).neighbor_up_right = PHTelem(e).neighbor_right(4);
        PHTelem(lastElem+6).neighbor_down_right = PHTelem(e).neighbor_right(1);
        PHTelem(lastElem+8).neighbor_down_right = PHTelem(e).neighbor_right(2);
        PHTelem(lastElem+2).neighbor_right_back = PHTelem(e).neighbor_right(2);
        PHTelem(lastElem+4).neighbor_right_front = PHTelem(e).neighbor_right(1);
        PHTelem(lastElem+6).neighbor_right_back = PHTelem(e).neighbor_right(4);
        PHTelem(lastElem+8).neighbor_right_front = PHTelem(e).neighbor_right(3);
        
    else
        PHTelem(lastElem+2).neighbor_right = PHTelem(e).neighbor_right;
        PHTelem(lastElem+4).neighbor_right = PHTelem(e).neighbor_right;
        PHTelem(lastElem+6).neighbor_right = PHTelem(e).neighbor_right;
        PHTelem(lastElem+8).neighbor_right = PHTelem(e).neighbor_right;
    end
    
    if (length(PHTelem(e).neighbor_up)==4)
        PHTelem(lastElem+5).neighbor_up = PHTelem(e).neighbor_up(1);
        PHTelem(lastElem+6).neighbor_up = PHTelem(e).neighbor_up(2);
        PHTelem(lastElem+7).neighbor_up = PHTelem(e).neighbor_up(3);
        PHTelem(lastElem+8).neighbor_up = PHTelem(e).neighbor_up(4);
        PHTelem(lastElem+5).neighbor_up_right = PHTelem(e).neighbor_up(2);
        PHTelem(lastElem+5).neighbor_up_back = PHTelem(e).neighbor_up(3);
        PHTelem(lastElem+6).neighbor_up_left = PHTelem(e).neighbor_up(1);
        PHTelem(lastElem+6).neighbor_up_back = PHTelem(e).neighbor_up(4);
        PHTelem(lastElem+7).neighbor_up_right = PHTelem(e).neighbor_up(4);
        PHTelem(lastElem+7).neighbor_up_front = PHTelem(e).neighbor_up(1);
        PHTelem(lastElem+8).neighbor_up_left = PHTelem(e).neighbor_up(3);
        PHTelem(lastElem+8).neighbor_up_front = PHTelem(e).neighbor_up(2);
    else
        PHTelem(lastElem+5).neighbor_up = PHTelem(e).neighbor_up;
        PHTelem(lastElem+6).neighbor_up = PHTelem(e).neighbor_up;
        PHTelem(lastElem+7).neighbor_up = PHTelem(e).neighbor_up;
        PHTelem(lastElem+8).neighbor_up = PHTelem(e).neighbor_up;
    end
    
    if (length(PHTelem(e).neighbor_left)==4)
        PHTelem(lastElem+1).neighbor_left = PHTelem(e).neighbor_left(1);
        PHTelem(lastElem+3).neighbor_left = PHTelem(e).neighbor_left(2);
        PHTelem(lastElem+5).neighbor_left = PHTelem(e).neighbor_left(3);
        PHTelem(lastElem+7).neighbor_left = PHTelem(e).neighbor_left(4);
        PHTelem(lastElem+1).neighbor_up_left = PHTelem(e).neighbor_left(3);
        PHTelem(lastElem+1).neighbor_left_back = PHTelem(e).neighbor_left(2);
        PHTelem(lastElem+3).neighbor_up_left = PHTelem(e).neighbor_left(4);
        PHTelem(lastElem+3).neighbor_left_front = PHTelem(e).neighbor_left(1);
        PHTelem(lastElem+5).neighbor_down_left = PHTelem(e).neighbor_left(1);
        PHTelem(lastElem+5).neighbor_left_back = PHTelem(e).neighbor_left(4);
        PHTelem(lastElem+7).neighbor_down_left = PHTelem(e).neighbor_left(2);
        PHTelem(lastElem+7).neighbor_left_front = PHTelem(e).neighbor_left(3);
    else
        PHTelem(lastElem+1).neighbor_left = PHTelem(e).neighbor_left;
        PHTelem(lastElem+3).neighbor_left = PHTelem(e).neighbor_left;
        PHTelem(lastElem+5).neighbor_left = PHTelem(e).neighbor_left;
        PHTelem(lastElem+7).neighbor_left = PHTelem(e).neighbor_left;
    end
    
    if (length(PHTelem(e).neighbor_front)==4)
        PHTelem(lastElem+1).neighbor_front = PHTelem(e).neighbor_front(1);
        PHTelem(lastElem+2).neighbor_front = PHTelem(e).neighbor_front(2);
        PHTelem(lastElem+5).neighbor_front = PHTelem(e).neighbor_front(3);
        PHTelem(lastElem+6).neighbor_front = PHTelem(e).neighbor_front(4);
        PHTelem(lastElem+1).neighbor_right_front = PHTelem(e).neighbor_front(2);
        PHTelem(lastElem+1).neighbor_up_front = PHTelem(e).neighbor_front(3);
        PHTelem(lastElem+2).neighbor_left_front = PHTelem(e).neighbor_front(1);
        PHTelem(lastElem+2).neighbor_up_front = PHTelem(e).neighbor_front(4);
        PHTelem(lastElem+5).neighbor_right_front = PHTelem(e).neighbor_front(4);
        PHTelem(lastElem+5).neighbor_down_front = PHTelem(e).neighbor_front(1);
        PHTelem(lastElem+6).neighbor_left_front = PHTelem(e).neighbor_front(3);
        PHTelem(lastElem+6).neighbor_down_front = PHTelem(e).neighbor_front(2);
    else
        PHTelem(lastElem+1).neighbor_front = PHTelem(e).neighbor_front;
        PHTelem(lastElem+2).neighbor_front = PHTelem(e).neighbor_front;
        PHTelem(lastElem+5).neighbor_front = PHTelem(e).neighbor_front;
        PHTelem(lastElem+6).neighbor_front = PHTelem(e).neighbor_front;
    end
    
    if (length(PHTelem(e).neighbor_back)==4)
        PHTelem(lastElem+3).neighbor_back = PHTelem(e).neighbor_back(1);
        PHTelem(lastElem+4).neighbor_back = PHTelem(e).neighbor_back(2);
        PHTelem(lastElem+7).neighbor_back = PHTelem(e).neighbor_back(3);
        PHTelem(lastElem+8).neighbor_back = PHTelem(e).neighbor_back(4);
        PHTelem(lastElem+3).neighbor_right_back = PHTelem(e).neighbor_back(2);
        PHTelem(lastElem+3).neighbor_up_back = PHTelem(e).neighbor_back(3);
        PHTelem(lastElem+4).neighbor_left_back = PHTelem(e).neighbor_back(1);
        PHTelem(lastElem+4).neighbor_up_back = PHTelem(e).neighbor_back(4);
        PHTelem(lastElem+7).neighbor_right_back = PHTelem(e).neighbor_back(4);
        PHTelem(lastElem+7).neighbor_down_back = PHTelem(e).neighbor_back(1);
        PHTelem(lastElem+8).neighbor_left_back = PHTelem(e).neighbor_back(3);
        PHTelem(lastElem+8).neighbor_down_back = PHTelem(e).neighbor_back(2);
    else
        PHTelem(lastElem+3).neighbor_back = PHTelem(e).neighbor_back;
        PHTelem(lastElem+4).neighbor_back = PHTelem(e).neighbor_back;
        PHTelem(lastElem+7).neighbor_back = PHTelem(e).neighbor_back;
        PHTelem(lastElem+8).neighbor_back = PHTelem(e).neighbor_back;
    end
    
    if (length(PHTelem(e).neighbor_up_left)==2)
        PHTelem(lastElem+5).neighbor_up_left = PHTelem(e).neighbor_up_left(1);
        PHTelem(lastElem+7).neighbor_up_left = PHTelem(e).neighbor_up_left(2);
    else
        PHTelem(lastElem+5).neighbor_up_left = PHTelem(e).neighbor_up_left;
        PHTelem(lastElem+7).neighbor_up_left = PHTelem(e).neighbor_up_left;
    end
    
    if (length(PHTelem(e).neighbor_down_left)==2)
        PHTelem(lastElem+1).neighbor_down_left = PHTelem(e).neighbor_down_left(1);
        PHTelem(lastElem+3).neighbor_down_left = PHTelem(e).neighbor_down_left(2);
    else
        PHTelem(lastElem+1).neighbor_down_left = PHTelem(e).neighbor_down_left;
        PHTelem(lastElem+3).neighbor_down_left = PHTelem(e).neighbor_down_left;
    end
    
    if (length(PHTelem(e).neighbor_up_right)==2)
        PHTelem(lastElem+6).neighbor_up_right = PHTelem(e).neighbor_up_right(1);
        PHTelem(lastElem+8).neighbor_up_right = PHTelem(e).neighbor_up_right(2);
    else
        PHTelem(lastElem+6).neighbor_up_right = PHTelem(e).neighbor_up_right;
        PHTelem(lastElem+8).neighbor_up_right = PHTelem(e).neighbor_up_right;
    end
    
    if (length(PHTelem(e).neighbor_down_right)==2)
        PHTelem(lastElem+2).neighbor_down_right = PHTelem(e).neighbor_down_right(1);
        PHTelem(lastElem+4).neighbor_down_right = PHTelem(e).neighbor_down_right(2);
    else
        PHTelem(lastElem+2).neighbor_down_right = PHTelem(e).neighbor_down_right;
        PHTelem(lastElem+4).neighbor_down_right = PHTelem(e).neighbor_down_right;
    end
    
    if (length(PHTelem(e).neighbor_up_front)==2)
        PHTelem(lastElem+5).neighbor_up_front = PHTelem(e).neighbor_up_front(1);
        PHTelem(lastElem+6).neighbor_up_front = PHTelem(e).neighbor_up_front(2);
    else
        PHTelem(lastElem+5).neighbor_up_front = PHTelem(e).neighbor_up_front;
        PHTelem(lastElem+6).neighbor_up_front = PHTelem(e).neighbor_up_front;
    end
    
    if (length(PHTelem(e).neighbor_down_front)==2)
        PHTelem(lastElem+1).neighbor_down_front = PHTelem(e).neighbor_down_front(1);
        PHTelem(lastElem+2).neighbor_down_front = PHTelem(e).neighbor_down_front(2);
    else
        PHTelem(lastElem+1).neighbor_down_front = PHTelem(e).neighbor_down_front;
        PHTelem(lastElem+2).neighbor_down_front = PHTelem(e).neighbor_down_front;
    end
    
    if (length(PHTelem(e).neighbor_up_back)==2)
        PHTelem(lastElem+7).neighbor_up_back = PHTelem(e).neighbor_up_back(1);
        PHTelem(lastElem+8).neighbor_up_back = PHTelem(e).neighbor_up_back(2);
    else
        PHTelem(lastElem+7).neighbor_up_back = PHTelem(e).neighbor_up_back;
        PHTelem(lastElem+8).neighbor_up_back = PHTelem(e).neighbor_up_back;
    end
    
    if (length(PHTelem(e).neighbor_down_back)==2)
        PHTelem(lastElem+3).neighbor_down_back = PHTelem(e).neighbor_down_back(1);
        PHTelem(lastElem+4).neighbor_down_back = PHTelem(e).neighbor_down_back(2);
    else
        PHTelem(lastElem+3).neighbor_down_back = PHTelem(e).neighbor_down_back;
        PHTelem(lastElem+4).neighbor_down_back = PHTelem(e).neighbor_down_back;
    end
    
    if (length(PHTelem(e).neighbor_left_front)==2)
        PHTelem(lastElem+1).neighbor_left_front = PHTelem(e).neighbor_left_front(1);
        PHTelem(lastElem+5).neighbor_left_front = PHTelem(e).neighbor_left_front(2);
    else
        PHTelem(lastElem+1).neighbor_left_front = PHTelem(e).neighbor_left_front;
        PHTelem(lastElem+5).neighbor_left_front = PHTelem(e).neighbor_left_front;
    end
    
    if (length(PHTelem(e).neighbor_right_front)==2)
        PHTelem(lastElem+2).neighbor_right_front = PHTelem(e).neighbor_right_front(1);
        PHTelem(lastElem+6).neighbor_right_front = PHTelem(e).neighbor_right_front(2);
    else
        PHTelem(lastElem+2).neighbor_right_front = PHTelem(e).neighbor_right_front;
        PHTelem(lastElem+6).neighbor_right_front = PHTelem(e).neighbor_right_front;
    end
    
    if (length(PHTelem(e).neighbor_left_back)==2)
        PHTelem(lastElem+3).neighbor_left_back = PHTelem(e).neighbor_left_back(1);
        PHTelem(lastElem+7).neighbor_left_back = PHTelem(e).neighbor_left_back(2);
    else
        PHTelem(lastElem+3).neighbor_left_back = PHTelem(e).neighbor_left_back;
        PHTelem(lastElem+7).neighbor_left_back = PHTelem(e).neighbor_left_back;
    end
    
    if (length(PHTelem(e).neighbor_right_back)==2)
        PHTelem(lastElem+4).neighbor_right_back = PHTelem(e).neighbor_right_back(1);
        PHTelem(lastElem+8).neighbor_right_back = PHTelem(e).neighbor_right_back(2);
    else
        PHTelem(lastElem+4).neighbor_right_back = PHTelem(e).neighbor_right_back;
        PHTelem(lastElem+8).neighbor_right_back = PHTelem(e).neighbor_right_back;
    end
    
    %handle the removed T-junctions
    dimBasisTemp = dimBasis;
    for ijunct = RTjunct
        parent_neighbor = zeros(1,3);
        %update the neighbors of the neighbor with children elements so we
        %get correct knotUln, knotUrn, knotVdn, knotVun      
        
        
        switch ijunct
            case 1
                front_neighbor = PHTelem(e).neighbor_front(1);
                parent_neighbor(1) = PHTelem(front_neighbor).parent;
                PHTelem(parent_neighbor(1)).neighbor_back = [lastElem+1, lastElem+2, lastElem+5, lastElem+6];
            case 2
                right_neighbor = PHTelem(e).neighbor_right(1);
                parent_neighbor(1) = PHTelem(right_neighbor).parent;
                PHTelem(parent_neighbor(1)).neighbor_left = [lastElem+2, lastElem+4, lastElem+6, lastElem+8];
            case 3
                back_neighbor = PHTelem(e).neighbor_back(1);
                parent_neighbor(1) = PHTelem(back_neighbor).parent;
                PHTelem(parent_neighbor(1)).neighbor_front = [lastElem+3, lastElem+4, lastElem+7, lastElem+8];
            case 4
                left_neighbor = PHTelem(e).neighbor_left(1);
                parent_neighbor(1) = PHTelem(left_neighbor).parent;
                PHTelem(parent_neighbor(1)).neighbor_right = [lastElem+1, lastElem+3, lastElem+5, lastElem+7];
            case 5
                down_neighbor = PHTelem(e).neighbor_down(1);
                parent_neighbor(1) = PHTelem(down_neighbor).parent;
                PHTelem(parent_neighbor(1)).neighbor_up = [lastElem+1, lastElem+2, lastElem+3, lastElem+4];
            case 6
                up_neighbor = PHTelem(e).neighbor_up(1);
                parent_neighbor(1) = PHTelem(up_neighbor).parent;
                PHTelem(parent_neighbor(1)).neighbor_down = [lastElem+5, lastElem+6, lastElem+7, lastElem+8];
            case 8
                if ~isempty(PHTelem(e).neighbor_up_left)
                    up_left_neighbor = PHTelem(e).neighbor_up_left(1);
                    parent_neighbor(1) = PHTelem(up_left_neighbor).parent;
                    PHTelem(parent_neighbor(1)).neighbor_down_right = [lastElem+5, lastElem+7];
                    PHTelem(PHTelem(PHTelem(e).neighbor_up(1)).parent).neighbor_down_left = PHTelem(e).neighbor_left([3,4]);
                    PHTelem(PHTelem(PHTelem(e).neighbor_left(1)).parent).neighbor_up_right = PHTelem(e).neighbor_up([1,3]);
                end
                if ~isempty(PHTelem(e).neighbor_up)
                    up_neighbor = PHTelem(e).neighbor_up(1);
                    parent_neighbor(2) = PHTelem(up_neighbor).parent;
                end
                if ~isempty(PHTelem(e).neighbor_left)
                    left_neighbor = PHTelem(e).neighbor_left(1);
                    parent_neighbor(3) = PHTelem(left_neighbor).parent;
                end
            case 9
                if ~isempty(PHTelem(e).neighbor_down_left)
                    down_left_neighbor = PHTelem(e).neighbor_down_left(1);
                    parent_neighbor(1) = PHTelem(down_left_neighbor).parent;
                    PHTelem(parent_neighbor(1)).neighbor_up_right = [lastElem+1, lastElem+3];
                    PHTelem(PHTelem(PHTelem(e).neighbor_down(1)).parent).neighbor_up_left = PHTelem(e).neighbor_left([1,2]);
                    PHTelem(PHTelem(PHTelem(e).neighbor_left(1)).parent).neighbor_down_right = PHTelem(e).neighbor_down([1,3]);
                end
                if ~isempty(PHTelem(e).neighbor_down)
                    down_neighbor = PHTelem(e).neighbor_down(1);
                    parent_neighbor(2) = PHTelem(down_neighbor).parent;
                end
                if ~isempty(PHTelem(e).neighbor_left)
                    left_neighbor = PHTelem(e).neighbor_left(1);
                    parent_neighbor(3) = PHTelem(left_neighbor).parent;
                end
            case 10
                if ~isempty(PHTelem(e).neighbor_up_right)
                    up_right_neighbor = PHTelem(e).neighbor_up_right(1);
                    parent_neighbor(1) = PHTelem(up_right_neighbor).parent;
                    PHTelem(parent_neighbor(1)).neighbor_down_left = [lastElem+6, lastElem+8];
                    PHTelem(PHTelem(PHTelem(e).neighbor_up(1)).parent).neighbor_down_right = PHTelem(e).neighbor_right([3,4]);
                    PHTelem(PHTelem(PHTelem(e).neighbor_right(1)).parent).neighbor_up_left = PHTelem(e).neighbor_up([2,4]);
                end
                if ~isempty(PHTelem(e).neighbor_up)
                    up_neighbor = PHTelem(e).neighbor_up(1);
                    parent_neighbor(2) = PHTelem(up_neighbor).parent;
                end
                if ~isempty(PHTelem(e).neighbor_right)
                    right_neighbor = PHTelem(e).neighbor_right(1);
                    parent_neighbor(3) = PHTelem(right_neighbor).parent;
                end
            case 11
                if ~isempty(PHTelem(e).neighbor_down_right)
                    down_right_neighbor = PHTelem(e).neighbor_down_right(1);
                    parent_neighbor(1) = PHTelem(down_right_neighbor).parent;
                    PHTelem(parent_neighbor(1)).neighbor_up_left = [lastElem+2, lastElem+4];
                    PHTelem(PHTelem(PHTelem(e).neighbor_down(1)).parent).neighbor_up_right = PHTelem(e).neighbor_right([1,3]);
                    PHTelem(PHTelem(PHTelem(e).neighbor_right(1)).parent).neighbor_down_left = PHTelem(e).neighbor_down([2,4]);
                end
                if ~isempty(PHTelem(e).neighbor_down)
                    down_neighbor = PHTelem(e).neighbor_down(1);
                    parent_neighbor(2) = PHTelem(down_neighbor).parent;
                end
                if ~isempty(PHTelem(e).neighbor_right)
                    right_neighbor = PHTelem(e).neighbor_right(1);
                    parent_neighbor(3) = PHTelem(right_neighbor).parent;
                end
            case 12
                if ~isempty(PHTelem(e).neighbor_up_front)
                    up_front_neighbor = PHTelem(e).neighbor_up_front(1);
                    parent_neighbor(1) = PHTelem(up_front_neighbor).parent;
                    PHTelem(parent_neighbor(1)).neighbor_down_back = [lastElem+5, lastElem+6];
                    PHTelem(PHTelem(PHTelem(e).neighbor_up(1)).parent).neighbor_down_front = PHTelem(e).neighbor_front([3,4]);
                    PHTelem(PHTelem(PHTelem(e).neighbor_front(1)).parent).neighbor_up_back = PHTelem(e).neighbor_up([1,2]);
                end
                if ~isempty(PHTelem(e).neighbor_up)
                    up_neighbor = PHTelem(e).neighbor_up(1);
                    parent_neighbor(2) = PHTelem(up_neighbor).parent;
                end
                if ~isempty(PHTelem(e).neighbor_front)
                    front_neighbor = PHTelem(e).neighbor_front(1);
                    parent_neighbor(3) = PHTelem(front_neighbor).parent;
                end
                
            case 13
                if ~isempty(PHTelem(e).neighbor_down_front)
                    down_front_neighbor = PHTelem(e).neighbor_down_front(1);
                    parent_neighbor(1) = PHTelem(down_front_neighbor).parent;
                    PHTelem(parent_neighbor(1)).neighbor_up_back = [lastElem+1, lastElem+2];
                    PHTelem(PHTelem(PHTelem(e).neighbor_down(1)).parent).neighbor_up_front = PHTelem(e).neighbor_front([1,2]);
                    PHTelem(PHTelem(PHTelem(e).neighbor_front(1)).parent).neighbor_down_back = PHTelem(e).neighbor_down([1,2]);
                end
                if ~isempty(PHTelem(e).neighbor_down)
                    down_neighbor = PHTelem(e).neighbor_down(1);
                    parent_neighbor(2) = PHTelem(down_neighbor).parent;
                end
                if ~isempty(PHTelem(e).neighbor_front)
                    front_neighbor = PHTelem(e).neighbor_front(1);
                    parent_neighbor(3) = PHTelem(front_neighbor).parent;
                end
            case 14
                if ~isempty(PHTelem(e).neighbor_up_back)
                    up_back_neighbor = PHTelem(e).neighbor_up_back(1);
                    parent_neighbor(1) = PHTelem(up_back_neighbor).parent;
                    PHTelem(parent_neighbor(1)).neighbor_down_front = [lastElem+7, lastElem+8];
                    PHTelem(PHTelem(PHTelem(e).neighbor_up(1)).parent).neighbor_down_back = PHTelem(e).neighbor_back([3,4]);
                    PHTelem(PHTelem(PHTelem(e).neighbor_back(1)).parent).neighbor_up_front = PHTelem(e).neighbor_up([3,4]);
                end
                if ~isempty(PHTelem(e).neighbor_up)
                    up_neighbor = PHTelem(e).neighbor_up(1);
                    parent_neighbor(2) = PHTelem(up_neighbor).parent;
                end
                if ~isempty(PHTelem(e).neighbor_back)
                    back_neighbor = PHTelem(e).neighbor_back(1);
                    parent_neighbor(3) = PHTelem(back_neighbor).parent;
                end
            case 15
                if ~isempty(PHTelem(e).neighbor_down_back)
                    down_back_neighbor = PHTelem(e).neighbor_down_back(1);
                    parent_neighbor(1) = PHTelem(down_back_neighbor).parent;
                    PHTelem(parent_neighbor(1)).neighbor_up_front = [lastElem+3, lastElem+4];
                    PHTelem(PHTelem(PHTelem(e).neighbor_down(1)).parent).neighbor_up_back = PHTelem(e).neighbor_back([1,2]);
                    PHTelem(PHTelem(PHTelem(e).neighbor_back(1)).parent).neighbor_down_front = PHTelem(e).neighbor_down([3,4]);
                end
                if ~isempty(PHTelem(e).neighbor_down)
                    down_neighbor = PHTelem(e).neighbor_down(1);
                    parent_neighbor(2) = PHTelem(down_neighbor).parent;
                end
                if ~isempty(PHTelem(e).neighbor_back)
                    back_neighbor = PHTelem(e).neighbor_back(1);
                    parent_neighbor(3) = PHTelem(back_neighbor).parent;
                end
            case 16
                if ~isempty(PHTelem(e).neighbor_left_front)
                    left_front_neighbor = PHTelem(e).neighbor_left_front(1);
                    parent_neighbor(1) = PHTelem(left_front_neighbor).parent;
                    PHTelem(parent_neighbor(1)).neighbor_right_back = [lastElem+1, lastElem+5];
                    PHTelem(PHTelem(PHTelem(e).neighbor_left(1)).parent).neighbor_right_front = PHTelem(e).neighbor_front([1,3]);
                    PHTelem(PHTelem(PHTelem(e).neighbor_front(1)).parent).neighbor_left_back = PHTelem(e).neighbor_left([1,3]);
                end
                if ~isempty(PHTelem(e).neighbor_left)
                    left_neighbor = PHTelem(e).neighbor_left(1);
                    parent_neighbor(2) = PHTelem(left_neighbor).parent;
                end
                if ~isempty(PHTelem(e).neighbor_front)
                    front_neighbor = PHTelem(e).neighbor_front(1);
                    parent_neighbor(3) = PHTelem(front_neighbor).parent;
                end
            case 17
                if ~isempty(PHTelem(e).neighbor_right_front)
                    right_front_neighbor = PHTelem(e).neighbor_right_front(1);
                    parent_neighbor(1) = PHTelem(right_front_neighbor).parent;
                    PHTelem(parent_neighbor(1)).neighbor_left_back = [lastElem+2, lastElem+6];
                    PHTelem(PHTelem(PHTelem(e).neighbor_right(1)).parent).neighbor_left_front = PHTelem(e).neighbor_front([2,4]);
                    PHTelem(PHTelem(PHTelem(e).neighbor_front(1)).parent).neighbor_right_back = PHTelem(e).neighbor_right([1,3]);
                end
                if ~isempty(PHTelem(e).neighbor_right)
                    right_neighbor = PHTelem(e).neighbor_right(1);
                    parent_neighbor(2) = PHTelem(right_neighbor).parent;
                end
                if ~isempty(PHTelem(e).neighbor_front)
                    front_neighbor = PHTelem(e).neighbor_front(1);
                    parent_neighbor(3) = PHTelem(front_neighbor).parent;
                end
            case 18
                if ~isempty(PHTelem(e).neighbor_left_back)
                    left_back_neighbor = PHTelem(e).neighbor_left_back(1);
                    parent_neighbor(1) = PHTelem(left_back_neighbor).parent;
                    PHTelem(parent_neighbor(1)).neighbor_right_front = [lastElem+3, lastElem+7];
                    PHTelem(PHTelem(PHTelem(e).neighbor_left(1)).parent).neighbor_right_back = PHTelem(e).neighbor_back([1,3]);
                    PHTelem(PHTelem(PHTelem(e).neighbor_back(1)).parent).neighbor_left_front = PHTelem(e).neighbor_left([2,4]);
                end
                if ~isempty(PHTelem(e).neighbor_left)
                    left_neighbor = PHTelem(e).neighbor_left(1);
                    parent_neighbor(2) = PHTelem(left_neighbor).parent;
                end
                if ~isempty(PHTelem(e).neighbor_back)
                    back_neighbor = PHTelem(e).neighbor_back(1);
                    parent_neighbor(3) = PHTelem(back_neighbor).parent;
                end
            case 19
                if ~isempty(PHTelem(e).neighbor_right_back)
                    right_back_neighbor = PHTelem(e).neighbor_right_back(1);
                    parent_neighbor(1) = PHTelem(right_back_neighbor).parent;
                    PHTelem(parent_neighbor(1)).neighbor_left_front = [lastElem+4, lastElem+8];
                    PHTelem(PHTelem(PHTelem(e).neighbor_right(1)).parent(1)).neighbor_left_back = PHTelem(e).neighbor_back([2,4]);
                    PHTelem(PHTelem(PHTelem(e).neighbor_back(1)).parent(1)).neighbor_right_front = PHTelem(e).neighbor_right([2,4]);
                    
                end
                if ~isempty(PHTelem(e).neighbor_right)
                    right_neighbor = PHTelem(e).neighbor_right(1);
                    parent_neighbor(2) = PHTelem(right_neighbor).parent;
                end
                if ~isempty(PHTelem(e).neighbor_back)
                    back_neighbor = PHTelem(e).neighbor_back(1);
                    parent_neighbor(3) = PHTelem(back_neighbor).parent;
                end
        end
        for parentIndex = 1:3
            
            if parent_neighbor(parentIndex) ~= 0
                
                [newBasisVertn,~,RTjunctn,knotUln, knotUrn, knotVdn, knotVun, knotWbn, knotWtn ] = checkNeighbors3D( PHTelem, parent_neighbor(parentIndex) );
                %update the Bezier extraction operators and element_nod indices
                C_n = PHTelem(parent_neighbor(parentIndex)).C;
                elmIn_n = PHTelem(parent_neighbor(parentIndex)).nodes;
                knotU1n = PHTelem(parent_neighbor(parentIndex)).vertex(1);
                knotU2n = PHTelem(parent_neighbor(parentIndex)).vertex(4);
                knotV1n = PHTelem(parent_neighbor(parentIndex)).vertex(2);
                knotV2n = PHTelem(parent_neighbor(parentIndex)).vertex(5);
                knotW1n = PHTelem(parent_neighbor(parentIndex)).vertex(3);
                knotW2n = PHTelem(parent_neighbor(parentIndex)).vertex(6);
                
                [bezt] = deCasteljau3dai( C_n, knotU1n, knotU2n, knotV1n, knotV2n, knotW1n, knotW2n, p, q, r, [RTjunctn, newBasisVertn, 7], knotUln, knotUrn, knotVdn, knotVun, knotWbn, knotWtn, elmIn_n, dimBasisTemp );
                
                switch ijunct
                    case 1
                        [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 3, elmIn_n, p, q, r, dimBasisTemp, reg );
                        front_neighbors = PHTelem(e).neighbor_front;
                        PHTelem(front_neighbors(1)).C = bezt.Ce3;
                        PHTelem(front_neighbors(2)).C = bezt.Ce4;
                        PHTelem(front_neighbors(3)).C = bezt.Ce7;
                        PHTelem(front_neighbors(4)).C = bezt.Ce8;
                        PHTelem(front_neighbors(1)).nodes(reg.north_mid) = newBasisSet{3}(reg.north_mid);
                        PHTelem(front_neighbors(1)).nodes(reg.ne_mid) = newBasisSet{3}(reg.ne_mid);
                        PHTelem(front_neighbors(1)).nodes(reg.north_top) = newBasisSet{3}(reg.north_top);
                        PHTelem(front_neighbors(1)).nodes(reg.ne_top) = newBasisSet{3}(reg.ne_top);
                        PHTelem(front_neighbors(2)).nodes(reg.north_mid) = newBasisSet{4}(reg.north_mid);
                        PHTelem(front_neighbors(2)).nodes(reg.nw_mid) = newBasisSet{4}(reg.nw_mid);
                        PHTelem(front_neighbors(2)).nodes(reg.north_top) = newBasisSet{4}(reg.north_top);
                        PHTelem(front_neighbors(2)).nodes(reg.nw_top) = newBasisSet{4}(reg.nw_top);
                        PHTelem(front_neighbors(3)).nodes(reg.north_low) = newBasisSet{7}(reg.north_low);
                        PHTelem(front_neighbors(3)).nodes(reg.ne_low) = newBasisSet{7}(reg.ne_low);
                        PHTelem(front_neighbors(3)).nodes(reg.north_mid) = newBasisSet{7}(reg.north_mid);
                        PHTelem(front_neighbors(3)).nodes(reg.ne_mid) = newBasisSet{7}(reg.ne_mid);
                        PHTelem(front_neighbors(4)).nodes(reg.north_low) = newBasisSet{8}(reg.north_low);
                        PHTelem(front_neighbors(4)).nodes(reg.nw_low) = newBasisSet{8}(reg.nw_low);
                        PHTelem(front_neighbors(4)).nodes(reg.north_mid) = newBasisSet{8}(reg.north_mid);
                        PHTelem(front_neighbors(4)).nodes(reg.nw_mid) = newBasisSet{8}(reg.nw_mid);
                        
                    case 2
                        [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 4, elmIn_n, p, q, r, dimBasisTemp, reg );
                        right_neighbors = PHTelem(e).neighbor_right;
                        PHTelem(right_neighbors(1)).C = bezt.Ce1;
                        PHTelem(right_neighbors(2)).C = bezt.Ce3;
                        PHTelem(right_neighbors(3)).C = bezt.Ce5;
                        PHTelem(right_neighbors(4)).C = bezt.Ce7;
                        PHTelem(right_neighbors(1)).nodes(reg.west_mid) = newBasisSet{1}(reg.west_mid);
                        PHTelem(right_neighbors(1)).nodes(reg.nw_mid) = newBasisSet{1}(reg.nw_mid);
                        PHTelem(right_neighbors(1)).nodes(reg.west_top) = newBasisSet{1}(reg.west_top);
                        PHTelem(right_neighbors(1)).nodes(reg.nw_top) = newBasisSet{1}(reg.nw_top);
                        PHTelem(right_neighbors(2)).nodes(reg.west_mid) = newBasisSet{3}(reg.west_mid);
                        PHTelem(right_neighbors(2)).nodes(reg.sw_mid) = newBasisSet{3}(reg.sw_mid);
                        PHTelem(right_neighbors(2)).nodes(reg.west_top) = newBasisSet{3}(reg.west_top);
                        PHTelem(right_neighbors(2)).nodes(reg.sw_top) = newBasisSet{3}(reg.sw_top);
                        PHTelem(right_neighbors(3)).nodes(reg.west_low) = newBasisSet{5}(reg.west_low);
                        PHTelem(right_neighbors(3)).nodes(reg.nw_low) = newBasisSet{5}(reg.nw_low);
                        PHTelem(right_neighbors(3)).nodes(reg.west_mid) = newBasisSet{5}(reg.west_mid);
                        PHTelem(right_neighbors(3)).nodes(reg.nw_mid) = newBasisSet{5}(reg.nw_mid);
                        PHTelem(right_neighbors(4)).nodes(reg.west_low) = newBasisSet{7}(reg.west_low);
                        PHTelem(right_neighbors(4)).nodes(reg.sw_low) = newBasisSet{7}(reg.sw_low);
                        PHTelem(right_neighbors(4)).nodes(reg.west_mid) = newBasisSet{7}(reg.west_mid);
                        PHTelem(right_neighbors(4)).nodes(reg.sw_mid) = newBasisSet{7}(reg.sw_mid);
                        
                    case 3
                        [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 1, elmIn_n, p, q, r, dimBasisTemp, reg );
                        back_neighbors = PHTelem(e).neighbor_back;
                        PHTelem(back_neighbors(1)).C = bezt.Ce1;
                        PHTelem(back_neighbors(2)).C = bezt.Ce2;
                        PHTelem(back_neighbors(3)).C = bezt.Ce5;
                        PHTelem(back_neighbors(4)).C = bezt.Ce6;
                        PHTelem(back_neighbors(1)).nodes(reg.south_mid) = newBasisSet{1}(reg.south_mid);
                        PHTelem(back_neighbors(1)).nodes(reg.se_mid) = newBasisSet{1}(reg.se_mid);
                        PHTelem(back_neighbors(1)).nodes(reg.south_top) = newBasisSet{1}(reg.south_top);
                        PHTelem(back_neighbors(1)).nodes(reg.se_top) = newBasisSet{1}(reg.se_top);
                        PHTelem(back_neighbors(2)).nodes(reg.south_mid) = newBasisSet{2}(reg.south_mid);
                        PHTelem(back_neighbors(2)).nodes(reg.sw_mid) = newBasisSet{2}(reg.sw_mid);
                        PHTelem(back_neighbors(2)).nodes(reg.south_top) = newBasisSet{2}(reg.south_top);
                        PHTelem(back_neighbors(2)).nodes(reg.sw_top) = newBasisSet{2}(reg.sw_top);
                        PHTelem(back_neighbors(3)).nodes(reg.south_low) = newBasisSet{5}(reg.south_low);
                        PHTelem(back_neighbors(3)).nodes(reg.se_low) = newBasisSet{5}(reg.se_low);
                        PHTelem(back_neighbors(3)).nodes(reg.south_mid) = newBasisSet{5}(reg.south_mid);
                        PHTelem(back_neighbors(3)).nodes(reg.se_mid) = newBasisSet{5}(reg.se_mid);
                        PHTelem(back_neighbors(4)).nodes(reg.south_low) = newBasisSet{6}(reg.south_low);
                        PHTelem(back_neighbors(4)).nodes(reg.sw_low) = newBasisSet{6}(reg.sw_low);
                        PHTelem(back_neighbors(4)).nodes(reg.south_mid) = newBasisSet{6}(reg.south_mid);
                        PHTelem(back_neighbors(4)).nodes(reg.sw_mid) = newBasisSet{6}(reg.sw_mid);
                        
                    case 4
                        [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 2, elmIn_n, p, q, r, dimBasisTemp, reg );
                        left_neighbors = PHTelem(e).neighbor_left;
                        PHTelem(left_neighbors(1)).C = bezt.Ce2;
                        PHTelem(left_neighbors(2)).C = bezt.Ce4;
                        PHTelem(left_neighbors(3)).C = bezt.Ce6;
                        PHTelem(left_neighbors(4)).C = bezt.Ce8;
                        PHTelem(left_neighbors(1)).nodes(reg.east_mid) = newBasisSet{2}(reg.east_mid);
                        PHTelem(left_neighbors(1)).nodes(reg.ne_mid) = newBasisSet{2}(reg.ne_mid);
                        PHTelem(left_neighbors(1)).nodes(reg.east_top) = newBasisSet{2}(reg.east_top);
                        PHTelem(left_neighbors(1)).nodes(reg.ne_top) = newBasisSet{2}(reg.ne_top);
                        PHTelem(left_neighbors(2)).nodes(reg.east_mid) = newBasisSet{4}(reg.east_mid);
                        PHTelem(left_neighbors(2)).nodes(reg.se_mid) = newBasisSet{4}(reg.se_mid);
                        PHTelem(left_neighbors(2)).nodes(reg.east_top) = newBasisSet{4}(reg.east_top);
                        PHTelem(left_neighbors(2)).nodes(reg.se_top) = newBasisSet{4}(reg.se_top);
                        PHTelem(left_neighbors(3)).nodes(reg.east_low) = newBasisSet{6}(reg.east_low);
                        PHTelem(left_neighbors(3)).nodes(reg.ne_low) = newBasisSet{6}(reg.ne_low);
                        PHTelem(left_neighbors(3)).nodes(reg.east_mid) = newBasisSet{6}(reg.east_mid);
                        PHTelem(left_neighbors(3)).nodes(reg.ne_mid) = newBasisSet{6}(reg.ne_mid);
                        PHTelem(left_neighbors(4)).nodes(reg.east_low) = newBasisSet{8}(reg.east_low);
                        PHTelem(left_neighbors(4)).nodes(reg.se_low) = newBasisSet{8}(reg.se_low);
                        PHTelem(left_neighbors(4)).nodes(reg.east_mid) = newBasisSet{8}(reg.east_mid);
                        PHTelem(left_neighbors(4)).nodes(reg.se_mid) = newBasisSet{8}(reg.se_mid);
                        
                    case 5
                        [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 6, elmIn_n, p, q, r, dimBasisTemp, reg );
                        down_neighbors = PHTelem(e).neighbor_down;
                        PHTelem(down_neighbors(1)).C = bezt.Ce5;
                        PHTelem(down_neighbors(2)).C = bezt.Ce6;
                        PHTelem(down_neighbors(3)).C = bezt.Ce7;
                        PHTelem(down_neighbors(4)).C = bezt.Ce8;
                        PHTelem(down_neighbors(1)).nodes(reg.center_top) = newBasisSet{5}(reg.center_top);
                        PHTelem(down_neighbors(1)).nodes(reg.east_top) = newBasisSet{5}(reg.east_top);
                        PHTelem(down_neighbors(1)).nodes(reg.north_top) = newBasisSet{5}(reg.north_top);
                        PHTelem(down_neighbors(1)).nodes(reg.ne_top) = newBasisSet{5}(reg.ne_top);
                        PHTelem(down_neighbors(2)).nodes(reg.center_top) = newBasisSet{6}(reg.center_top);
                        PHTelem(down_neighbors(2)).nodes(reg.west_top) = newBasisSet{6}(reg.west_top);
                        PHTelem(down_neighbors(2)).nodes(reg.north_top) = newBasisSet{6}(reg.north_top);
                        PHTelem(down_neighbors(2)).nodes(reg.nw_top) = newBasisSet{6}(reg.nw_top);
                        PHTelem(down_neighbors(3)).nodes(reg.center_top) = newBasisSet{7}(reg.center_top);
                        PHTelem(down_neighbors(3)).nodes(reg.south_top) = newBasisSet{7}(reg.south_top);
                        PHTelem(down_neighbors(3)).nodes(reg.east_top) = newBasisSet{7}(reg.east_top);
                        PHTelem(down_neighbors(3)).nodes(reg.se_top) = newBasisSet{7}(reg.se_top);
                        PHTelem(down_neighbors(4)).nodes(reg.center_top) = newBasisSet{8}(reg.center_top);
                        PHTelem(down_neighbors(4)).nodes(reg.west_top) = newBasisSet{8}(reg.west_top);
                        PHTelem(down_neighbors(4)).nodes(reg.south_top) = newBasisSet{8}(reg.south_top);
                        PHTelem(down_neighbors(4)).nodes(reg.sw_top) = newBasisSet{8}(reg.sw_top);
                        
                    case 6
                        [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 5, elmIn_n, p, q, r, dimBasisTemp, reg );
                        up_neighbors = PHTelem(e).neighbor_up;
                        PHTelem(up_neighbors(1)).C = bezt.Ce1;
                        PHTelem(up_neighbors(2)).C = bezt.Ce2;
                        PHTelem(up_neighbors(3)).C = bezt.Ce3;
                        PHTelem(up_neighbors(4)).C = bezt.Ce4;
                        PHTelem(up_neighbors(1)).nodes(reg.center_low) = newBasisSet{1}(reg.center_low);
                        PHTelem(up_neighbors(1)).nodes(reg.east_low) = newBasisSet{1}(reg.east_low);
                        PHTelem(up_neighbors(1)).nodes(reg.north_low) = newBasisSet{1}(reg.north_low);
                        PHTelem(up_neighbors(1)).nodes(reg.ne_low) = newBasisSet{1}(reg.ne_low);
                        PHTelem(up_neighbors(2)).nodes(reg.center_low) = newBasisSet{2}(reg.center_low);
                        PHTelem(up_neighbors(2)).nodes(reg.west_low) = newBasisSet{2}(reg.west_low);
                        PHTelem(up_neighbors(2)).nodes(reg.north_low) = newBasisSet{2}(reg.north_low);
                        PHTelem(up_neighbors(2)).nodes(reg.nw_low) = newBasisSet{2}(reg.nw_low);
                        PHTelem(up_neighbors(3)).nodes(reg.center_low) = newBasisSet{3}(reg.center_low);
                        PHTelem(up_neighbors(3)).nodes(reg.south_low) = newBasisSet{3}(reg.south_low);
                        PHTelem(up_neighbors(3)).nodes(reg.east_low) = newBasisSet{3}(reg.east_low);
                        PHTelem(up_neighbors(3)).nodes(reg.se_low) = newBasisSet{3}(reg.se_low);
                        PHTelem(up_neighbors(4)).nodes(reg.center_low) = newBasisSet{4}(reg.center_low);
                        PHTelem(up_neighbors(4)).nodes(reg.west_low) = newBasisSet{4}(reg.west_low);
                        PHTelem(up_neighbors(4)).nodes(reg.south_low) = newBasisSet{4}(reg.south_low);
                        PHTelem(up_neighbors(4)).nodes(reg.sw_low) = newBasisSet{4}(reg.sw_low);
                        
                    case 8
                        if ~isempty(PHTelem(e).neighbor_up_left) && (parentIndex==1)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 11, elmIn_n, p, q, r, dimBasisTemp, reg );
                            up_left_neighbors = PHTelem(e).neighbor_up_left;
                            PHTelem(up_left_neighbors(1)).C = bezt.Ce2;
                            PHTelem(up_left_neighbors(2)).C = bezt.Ce4;
                            PHTelem(up_left_neighbors(1)).nodes(reg.east_low) = newBasisSet{2}(reg.east_low);
                            PHTelem(up_left_neighbors(1)).nodes(reg.ne_low) = newBasisSet{2}(reg.ne_low);
                            PHTelem(up_left_neighbors(2)).nodes(reg.se_low) = newBasisSet{4}(reg.se_low);
                            PHTelem(up_left_neighbors(2)).nodes(reg.east_low) = newBasisSet{4}(reg.east_low);
                        end
                        if ~isempty(PHTelem(e).neighbor_up) && (parentIndex==2)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 9, elmIn_n, p, q, r, dimBasisTemp, reg );
                            up_left_neighbors = PHTelem(e).neighbor_up;
                            PHTelem(up_left_neighbors(1)).C = bezt.Ce1;
                            PHTelem(up_left_neighbors(3)).C = bezt.Ce3;
                            PHTelem(up_left_neighbors(1)).nodes(reg.west_low) = newBasisSet{1}(reg.west_low);
                            PHTelem(up_left_neighbors(1)).nodes(reg.nw_low) = newBasisSet{1}(reg.nw_low);
                            PHTelem(up_left_neighbors(3)).nodes(reg.sw_low) = newBasisSet{3}(reg.sw_low);
                            PHTelem(up_left_neighbors(3)).nodes(reg.west_low) = newBasisSet{3}(reg.west_low);
                        end
                        if ~isempty(PHTelem(e).neighbor_left) && (parentIndex==3)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 10, elmIn_n, p, q, r, dimBasisTemp, reg );
                            up_left_neighbors = PHTelem(e).neighbor_left;
                            PHTelem(up_left_neighbors(3)).C = bezt.Ce6;
                            PHTelem(up_left_neighbors(4)).C = bezt.Ce8;
                            PHTelem(up_left_neighbors(3)).nodes(reg.east_top) = newBasisSet{6}(reg.east_top);
                            PHTelem(up_left_neighbors(3)).nodes(reg.ne_top) = newBasisSet{6}(reg.ne_top);
                            PHTelem(up_left_neighbors(4)).nodes(reg.se_top) = newBasisSet{8}(reg.se_top);
                            PHTelem(up_left_neighbors(4)).nodes(reg.east_top) = newBasisSet{8}(reg.east_top);
                        end
                        
                    case 9
                        if ~isempty(PHTelem(e).neighbor_down_left) && (parentIndex==1)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 10, elmIn_n, p, q, r, dimBasisTemp, reg );
                            down_left_neighbors = PHTelem(e).neighbor_down_left;
                            PHTelem(down_left_neighbors(1)).C = bezt.Ce6;
                            PHTelem(down_left_neighbors(2)).C = bezt.Ce8;
                            PHTelem(down_left_neighbors(1)).nodes(reg.east_top) = newBasisSet{6}(reg.east_top);
                            PHTelem(down_left_neighbors(1)).nodes(reg.ne_top) = newBasisSet{6}(reg.ne_top);
                            PHTelem(down_left_neighbors(2)).nodes(reg.se_top) = newBasisSet{8}(reg.se_top);
                            PHTelem(down_left_neighbors(2)).nodes(reg.east_top) = newBasisSet{8}(reg.east_top);
                        end
                        if ~isempty(PHTelem(e).neighbor_down) && (parentIndex==2)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 8, elmIn_n, p, q, r, dimBasisTemp, reg );
                            down_left_neighbors = PHTelem(e).neighbor_down;
                            PHTelem(down_left_neighbors(1)).C = bezt.Ce5;
                            PHTelem(down_left_neighbors(3)).C = bezt.Ce7;
                            PHTelem(down_left_neighbors(1)).nodes(reg.west_top) = newBasisSet{5}(reg.west_top);
                            PHTelem(down_left_neighbors(1)).nodes(reg.nw_top) = newBasisSet{5}(reg.nw_top);
                            PHTelem(down_left_neighbors(3)).nodes(reg.sw_top) = newBasisSet{7}(reg.sw_top);
                            PHTelem(down_left_neighbors(3)).nodes(reg.west_top) = newBasisSet{7}(reg.west_top);
                        end
                        if ~isempty(PHTelem(e).neighbor_left) && (parentIndex==3)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 11, elmIn_n, p, q, r, dimBasisTemp, reg );
                            down_left_neighbors = PHTelem(e).neighbor_left;
                            PHTelem(down_left_neighbors(1)).C = bezt.Ce2;
                            PHTelem(down_left_neighbors(2)).C = bezt.Ce4;
                            PHTelem(down_left_neighbors(1)).nodes(reg.east_low) = newBasisSet{2}(reg.east_low);
                            PHTelem(down_left_neighbors(1)).nodes(reg.ne_low) = newBasisSet{2}(reg.ne_low);
                            PHTelem(down_left_neighbors(2)).nodes(reg.se_low) = newBasisSet{4}(reg.se_low);
                            PHTelem(down_left_neighbors(2)).nodes(reg.east_low) = newBasisSet{4}(reg.east_low);
                        end
                        
                    case 10
                        if ~isempty(PHTelem(e).neighbor_up_right) && (parentIndex==1)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 9, elmIn_n, p, q, r, dimBasisTemp, reg );
                            up_right_neighbors = PHTelem(e).neighbor_up_right;
                            PHTelem(up_right_neighbors(1)).C = bezt.Ce1;
                            PHTelem(up_right_neighbors(2)).C = bezt.Ce3;
                            PHTelem(up_right_neighbors(1)).nodes(reg.west_low) = newBasisSet{1}(reg.west_low);
                            PHTelem(up_right_neighbors(1)).nodes(reg.nw_low) = newBasisSet{1}(reg.nw_low);
                            PHTelem(up_right_neighbors(2)).nodes(reg.sw_low) = newBasisSet{3}(reg.sw_low);
                            PHTelem(up_right_neighbors(2)).nodes(reg.west_low) = newBasisSet{3}(reg.west_low);
                        end
                        if ~isempty(PHTelem(e).neighbor_up) && (parentIndex==2)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 11, elmIn_n, p, q, r, dimBasisTemp, reg );
                            up_right_neighbors = PHTelem(e).neighbor_up;
                            PHTelem(up_right_neighbors(2)).C = bezt.Ce2;
                            PHTelem(up_right_neighbors(4)).C = bezt.Ce4;
                            PHTelem(up_right_neighbors(2)).nodes(reg.east_low) = newBasisSet{2}(reg.east_low);
                            PHTelem(up_right_neighbors(2)).nodes(reg.ne_low) = newBasisSet{2}(reg.ne_low);
                            PHTelem(up_right_neighbors(4)).nodes(reg.se_low) = newBasisSet{4}(reg.se_low);
                            PHTelem(up_right_neighbors(4)).nodes(reg.east_low) = newBasisSet{4}(reg.east_low);
                        end
                        if ~isempty(PHTelem(e).neighbor_right) && (parentIndex==3)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 8, elmIn_n, p, q, r, dimBasisTemp, reg );
                            up_right_neighbors = PHTelem(e).neighbor_right;
                            PHTelem(up_right_neighbors(3)).C = bezt.Ce5;
                            PHTelem(up_right_neighbors(4)).C = bezt.Ce7;
                            PHTelem(up_right_neighbors(3)).nodes(reg.west_top) = newBasisSet{5}(reg.west_top);
                            PHTelem(up_right_neighbors(3)).nodes(reg.nw_top) = newBasisSet{5}(reg.nw_top);
                            PHTelem(up_right_neighbors(4)).nodes(reg.sw_top) = newBasisSet{7}(reg.sw_top);
                            PHTelem(up_right_neighbors(4)).nodes(reg.west_top) = newBasisSet{7}(reg.west_top);
                        end
                        
                    case 11
                        if ~isempty(PHTelem(e).neighbor_down_right) && (parentIndex==1)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 8, elmIn_n, p, q, r, dimBasisTemp, reg );
                            down_right_neighbors = PHTelem(e).neighbor_down_right;
                            PHTelem(down_right_neighbors(1)).C = bezt.Ce5;
                            PHTelem(down_right_neighbors(2)).C = bezt.Ce7;
                            PHTelem(down_right_neighbors(1)).nodes(reg.west_top) = newBasisSet{5}(reg.west_top);
                            PHTelem(down_right_neighbors(1)).nodes(reg.nw_top) = newBasisSet{5}(reg.nw_top);
                            PHTelem(down_right_neighbors(2)).nodes(reg.sw_top) = newBasisSet{7}(reg.sw_top);
                            PHTelem(down_right_neighbors(2)).nodes(reg.west_top) = newBasisSet{7}(reg.west_top);
                        end
                        if ~isempty(PHTelem(e).neighbor_down) && (parentIndex==2)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 10, elmIn_n, p, q, r, dimBasisTemp, reg );
                            down_right_neighbors = PHTelem(e).neighbor_down;
                            PHTelem(down_right_neighbors(2)).C = bezt.Ce6;
                            PHTelem(down_right_neighbors(4)).C = bezt.Ce8;
                            PHTelem(down_right_neighbors(2)).nodes(reg.east_top) = newBasisSet{6}(reg.east_top);
                            PHTelem(down_right_neighbors(2)).nodes(reg.ne_top) = newBasisSet{6}(reg.ne_top);
                            PHTelem(down_right_neighbors(4)).nodes(reg.se_top) = newBasisSet{8}(reg.se_top);
                            PHTelem(down_right_neighbors(4)).nodes(reg.east_top) = newBasisSet{8}(reg.east_top);
                        end
                        if ~isempty(PHTelem(e).neighbor_right) && (parentIndex==3)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 9, elmIn_n, p, q, r, dimBasisTemp, reg );
                            down_right_neighbors = PHTelem(e).neighbor_right;
                            PHTelem(down_right_neighbors(1)).C = bezt.Ce1;
                            PHTelem(down_right_neighbors(2)).C = bezt.Ce3;
                            PHTelem(down_right_neighbors(1)).nodes(reg.west_low) = newBasisSet{1}(reg.west_low);
                            PHTelem(down_right_neighbors(1)).nodes(reg.nw_low) = newBasisSet{1}(reg.nw_low);
                            PHTelem(down_right_neighbors(2)).nodes(reg.sw_low) = newBasisSet{3}(reg.sw_low);
                            PHTelem(down_right_neighbors(2)).nodes(reg.west_low) = newBasisSet{3}(reg.west_low);
                        end
                        
                    case 12
                        if ~isempty(PHTelem(e).neighbor_up_front) && (parentIndex==1)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 15, elmIn_n, p, q, r, dimBasisTemp, reg );
                            up_front_neighbors = PHTelem(e).neighbor_up_front;
                            PHTelem(up_front_neighbors(1)).C = bezt.Ce3;
                            PHTelem(up_front_neighbors(2)).C = bezt.Ce4;
                            PHTelem(up_front_neighbors(1)).nodes(reg.north_low) = newBasisSet{3}(reg.north_low);
                            PHTelem(up_front_neighbors(1)).nodes(reg.ne_low) = newBasisSet{3}(reg.ne_low);
                            PHTelem(up_front_neighbors(2)).nodes(reg.nw_low) = newBasisSet{4}(reg.nw_low);
                            PHTelem(up_front_neighbors(2)).nodes(reg.north_low) = newBasisSet{4}(reg.north_low);
                        end
                        if ~isempty(PHTelem(e).neighbor_up) && (parentIndex==2)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 13, elmIn_n, p, q, r, dimBasisTemp, reg );
                            up_front_neighbors = PHTelem(e).neighbor_up;
                            PHTelem(up_front_neighbors(1)).C = bezt.Ce1;
                            PHTelem(up_front_neighbors(2)).C = bezt.Ce2;
                            PHTelem(up_front_neighbors(1)).nodes(reg.south_low) = newBasisSet{1}(reg.south_low);
                            PHTelem(up_front_neighbors(1)).nodes(reg.se_low) = newBasisSet{1}(reg.se_low);
                            PHTelem(up_front_neighbors(2)).nodes(reg.sw_low) = newBasisSet{2}(reg.sw_low);
                            PHTelem(up_front_neighbors(2)).nodes(reg.south_low) = newBasisSet{2}(reg.south_low);
                        end
                        if ~isempty(PHTelem(e).neighbor_front) && (parentIndex==3)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 14, elmIn_n, p, q, r, dimBasisTemp, reg );
                            up_front_neighbors = PHTelem(e).neighbor_front;
                            PHTelem(up_front_neighbors(3)).C = bezt.Ce7;
                            PHTelem(up_front_neighbors(4)).C = bezt.Ce8;
                            PHTelem(up_front_neighbors(3)).nodes(reg.north_top) = newBasisSet{7}(reg.north_top);
                            PHTelem(up_front_neighbors(3)).nodes(reg.ne_top) = newBasisSet{7}(reg.ne_top);
                            PHTelem(up_front_neighbors(4)).nodes(reg.nw_top) = newBasisSet{8}(reg.nw_top);
                            PHTelem(up_front_neighbors(4)).nodes(reg.north_top) = newBasisSet{8}(reg.north_top);
                        end
                        
                    case 13
                        if ~isempty(PHTelem(e).neighbor_down_front) && (parentIndex==1)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 14, elmIn_n, p, q, r, dimBasisTemp, reg );
                            down_front_neighbors = PHTelem(e).neighbor_down_front;
                            PHTelem(down_front_neighbors(1)).C = bezt.Ce7;
                            PHTelem(down_front_neighbors(2)).C = bezt.Ce8;
                            PHTelem(down_front_neighbors(1)).nodes(reg.north_top) = newBasisSet{7}(reg.north_top);
                            PHTelem(down_front_neighbors(1)).nodes(reg.ne_top) = newBasisSet{7}(reg.ne_top);
                            PHTelem(down_front_neighbors(2)).nodes(reg.nw_top) = newBasisSet{8}(reg.nw_top);
                            PHTelem(down_front_neighbors(2)).nodes(reg.north_top) = newBasisSet{8}(reg.north_top);
                        end
                        
                        if ~isempty(PHTelem(e).neighbor_down) && (parentIndex==2)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 12, elmIn_n, p, q, r, dimBasisTemp, reg );
                            down_front_neighbors = PHTelem(e).neighbor_down;
                            PHTelem(down_front_neighbors(1)).C = bezt.Ce5;
                            PHTelem(down_front_neighbors(2)).C = bezt.Ce6;
                            PHTelem(down_front_neighbors(1)).nodes(reg.south_top) = newBasisSet{5}(reg.south_top);
                            PHTelem(down_front_neighbors(1)).nodes(reg.se_top) = newBasisSet{5}(reg.se_top);
                            PHTelem(down_front_neighbors(2)).nodes(reg.sw_top) = newBasisSet{6}(reg.sw_top);
                            PHTelem(down_front_neighbors(2)).nodes(reg.south_top) = newBasisSet{6}(reg.south_top);
                        end
                        
                        if ~isempty(PHTelem(e).neighbor_front) && (parentIndex==3)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 15, elmIn_n, p, q, r, dimBasisTemp, reg );
                            down_front_neighbors = PHTelem(e).neighbor_front;
                            PHTelem(down_front_neighbors(1)).C = bezt.Ce3;
                            PHTelem(down_front_neighbors(2)).C = bezt.Ce4;
                            PHTelem(down_front_neighbors(1)).nodes(reg.north_low) = newBasisSet{3}(reg.north_low);
                            PHTelem(down_front_neighbors(1)).nodes(reg.ne_low) = newBasisSet{3}(reg.ne_low);
                            PHTelem(down_front_neighbors(2)).nodes(reg.nw_low) = newBasisSet{4}(reg.nw_low);
                            PHTelem(down_front_neighbors(2)).nodes(reg.north_low) = newBasisSet{4}(reg.north_low);
                        end
                        
                        
                    case 14
                        if ~isempty(PHTelem(e).neighbor_up_back) && (parentIndex==1)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 13, elmIn_n, p, q, r, dimBasisTemp, reg );
                            up_back_neighbors = PHTelem(e).neighbor_up_back;
                            PHTelem(up_back_neighbors(1)).C = bezt.Ce1;
                            PHTelem(up_back_neighbors(2)).C = bezt.Ce2;
                            PHTelem(up_back_neighbors(1)).nodes(reg.south_low) = newBasisSet{1}(reg.south_low);
                            PHTelem(up_back_neighbors(1)).nodes(reg.se_low) = newBasisSet{1}(reg.se_low);
                            PHTelem(up_back_neighbors(2)).nodes(reg.sw_low) = newBasisSet{2}(reg.sw_low);
                            PHTelem(up_back_neighbors(2)).nodes(reg.south_low) = newBasisSet{2}(reg.south_low);
                        end
                        if ~isempty(PHTelem(e).neighbor_up) && (parentIndex==2)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 15, elmIn_n, p, q, r, dimBasisTemp, reg );
                            up_back_neighbors = PHTelem(e).neighbor_up;
                            PHTelem(up_back_neighbors(3)).C = bezt.Ce3;
                            PHTelem(up_back_neighbors(4)).C = bezt.Ce4;
                            PHTelem(up_back_neighbors(3)).nodes(reg.north_low) = newBasisSet{3}(reg.north_low);
                            PHTelem(up_back_neighbors(3)).nodes(reg.ne_low) = newBasisSet{3}(reg.ne_low);
                            PHTelem(up_back_neighbors(4)).nodes(reg.nw_low) = newBasisSet{4}(reg.nw_low);
                            PHTelem(up_back_neighbors(4)).nodes(reg.north_low) = newBasisSet{4}(reg.north_low);
                        end
                        if ~isempty(PHTelem(e).neighbor_back) && (parentIndex==3)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 12, elmIn_n, p, q, r, dimBasisTemp, reg );
                            up_back_neighbors = PHTelem(e).neighbor_back;
                            PHTelem(up_back_neighbors(3)).C = bezt.Ce5;
                            PHTelem(up_back_neighbors(4)).C = bezt.Ce6;
                            PHTelem(up_back_neighbors(3)).nodes(reg.south_top) = newBasisSet{5}(reg.south_top);
                            PHTelem(up_back_neighbors(3)).nodes(reg.se_top) = newBasisSet{5}(reg.se_top);
                            PHTelem(up_back_neighbors(4)).nodes(reg.sw_top) = newBasisSet{6}(reg.sw_top);
                            PHTelem(up_back_neighbors(4)).nodes(reg.south_top) = newBasisSet{6}(reg.south_top);
                        end
                        
                    case 15
                        if ~isempty(PHTelem(e).neighbor_down_back) && (parentIndex==1)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 12, elmIn_n, p, q, r, dimBasisTemp, reg );
                            down_back_neighbors = PHTelem(e).neighbor_down_back;
                            PHTelem(down_back_neighbors(1)).C = bezt.Ce5;
                            PHTelem(down_back_neighbors(2)).C = bezt.Ce6;
                            PHTelem(down_back_neighbors(1)).nodes(reg.south_top) = newBasisSet{5}(reg.south_top);
                            PHTelem(down_back_neighbors(1)).nodes(reg.se_top) = newBasisSet{5}(reg.se_top);
                            PHTelem(down_back_neighbors(2)).nodes(reg.sw_top) = newBasisSet{6}(reg.sw_top);
                            PHTelem(down_back_neighbors(2)).nodes(reg.south_top) = newBasisSet{6}(reg.south_top);
                        end
                        if ~isempty(PHTelem(e).neighbor_down) && (parentIndex==2)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 14, elmIn_n, p, q, r, dimBasisTemp, reg );
                            down_back_neighbors = PHTelem(e).neighbor_down;
                            PHTelem(down_back_neighbors(3)).C = bezt.Ce7;
                            PHTelem(down_back_neighbors(4)).C = bezt.Ce8;
                            PHTelem(down_back_neighbors(3)).nodes(reg.north_top) = newBasisSet{7}(reg.north_top);
                            PHTelem(down_back_neighbors(3)).nodes(reg.ne_top) = newBasisSet{7}(reg.ne_top);
                            PHTelem(down_back_neighbors(4)).nodes(reg.nw_top) = newBasisSet{8}(reg.nw_top);
                            PHTelem(down_back_neighbors(4)).nodes(reg.north_top) = newBasisSet{8}(reg.north_top);
                        end
                        if ~isempty(PHTelem(e).neighbor_back) && (parentIndex==3)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 13, elmIn_n, p, q, r, dimBasisTemp, reg );
                            down_back_neighbors = PHTelem(e).neighbor_back;
                            PHTelem(down_back_neighbors(1)).C = bezt.Ce1;
                            PHTelem(down_back_neighbors(2)).C = bezt.Ce2;
                            PHTelem(down_back_neighbors(1)).nodes(reg.south_low) = newBasisSet{1}(reg.south_low);
                            PHTelem(down_back_neighbors(1)).nodes(reg.se_low) = newBasisSet{1}(reg.se_low);
                            PHTelem(down_back_neighbors(2)).nodes(reg.sw_low) = newBasisSet{2}(reg.sw_low);
                            PHTelem(down_back_neighbors(2)).nodes(reg.south_low) = newBasisSet{2}(reg.south_low);
                        end
                        
                    case 16
                        if ~isempty(PHTelem(e).neighbor_left_front) && (parentIndex==1)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 19, elmIn_n, p, q, r, dimBasisTemp, reg );
                            left_front_neighbors = PHTelem(e).neighbor_left_front;
                            PHTelem(left_front_neighbors(1)).C = bezt.Ce4;
                            PHTelem(left_front_neighbors(2)).C = bezt.Ce8;
                            PHTelem(left_front_neighbors(1)).nodes(reg.ne_mid) = newBasisSet{4}(reg.ne_mid);
                            PHTelem(left_front_neighbors(1)).nodes(reg.ne_top) = newBasisSet{4}(reg.ne_top);
                            PHTelem(left_front_neighbors(2)).nodes(reg.ne_low) = newBasisSet{8}(reg.ne_low);
                            PHTelem(left_front_neighbors(2)).nodes(reg.ne_mid) = newBasisSet{8}(reg.ne_mid);
                        end
                        if ~isempty(PHTelem(e).neighbor_left) && (parentIndex==2)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 17, elmIn_n, p, q, r, dimBasisTemp, reg );
                            left_front_neighbors = PHTelem(e).neighbor_left;
                            PHTelem(left_front_neighbors(1)).C = bezt.Ce2;
                            PHTelem(left_front_neighbors(3)).C = bezt.Ce6;
                            PHTelem(left_front_neighbors(1)).nodes(reg.se_mid) = newBasisSet{2}(reg.se_mid);
                            PHTelem(left_front_neighbors(1)).nodes(reg.se_top) = newBasisSet{2}(reg.se_top);
                            PHTelem(left_front_neighbors(3)).nodes(reg.se_low) = newBasisSet{6}(reg.se_low);
                            PHTelem(left_front_neighbors(3)).nodes(reg.se_mid) = newBasisSet{6}(reg.se_mid);
                        end
                        if ~isempty(PHTelem(e).neighbor_front) && (parentIndex==3)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 18, elmIn_n, p, q, r, dimBasisTemp, reg );
                            left_front_neighbors = PHTelem(e).neighbor_front;
                            PHTelem(left_front_neighbors(1)).C = bezt.Ce3;
                            PHTelem(left_front_neighbors(3)).C = bezt.Ce7;
                            PHTelem(left_front_neighbors(1)).nodes(reg.nw_mid) = newBasisSet{3}(reg.nw_mid);
                            PHTelem(left_front_neighbors(1)).nodes(reg.nw_top) = newBasisSet{3}(reg.nw_top);
                            PHTelem(left_front_neighbors(3)).nodes(reg.nw_low) = newBasisSet{7}(reg.nw_low);
                            PHTelem(left_front_neighbors(3)).nodes(reg.nw_mid) = newBasisSet{7}(reg.nw_mid);
                        end
                        
                    case 17
                        if ~isempty(PHTelem(e).neighbor_right_front) && (parentIndex==1)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 18, elmIn_n, p, q, r, dimBasisTemp, reg );
                            right_front_neighbors = PHTelem(e).neighbor_right_front;
                            PHTelem(right_front_neighbors(1)).C = bezt.Ce3;
                            PHTelem(right_front_neighbors(2)).C = bezt.Ce7;
                            PHTelem(right_front_neighbors(1)).nodes(reg.nw_mid) = newBasisSet{3}(reg.nw_mid);
                            PHTelem(right_front_neighbors(1)).nodes(reg.nw_top) = newBasisSet{3}(reg.nw_top);
                            PHTelem(right_front_neighbors(2)).nodes(reg.nw_low) = newBasisSet{7}(reg.nw_low);
                            PHTelem(right_front_neighbors(2)).nodes(reg.nw_mid) = newBasisSet{7}(reg.nw_mid);
                        end
                        if ~isempty(PHTelem(e).neighbor_right) && (parentIndex==2)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 16, elmIn_n, p, q, r, dimBasisTemp, reg );
                            right_front_neighbors = PHTelem(e).neighbor_right;
                            PHTelem(right_front_neighbors(1)).C = bezt.Ce1;
                            PHTelem(right_front_neighbors(3)).C = bezt.Ce5;
                            PHTelem(right_front_neighbors(1)).nodes(reg.sw_mid) = newBasisSet{1}(reg.sw_mid);
                            PHTelem(right_front_neighbors(1)).nodes(reg.sw_top) = newBasisSet{1}(reg.sw_top);
                            PHTelem(right_front_neighbors(3)).nodes(reg.sw_low) = newBasisSet{5}(reg.sw_low);
                            PHTelem(right_front_neighbors(3)).nodes(reg.sw_mid) = newBasisSet{5}(reg.sw_mid);
                        end
                        if ~isempty(PHTelem(e).neighbor_front) && (parentIndex==3)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 19, elmIn_n, p, q, r, dimBasisTemp, reg );
                            right_front_neighbors = PHTelem(e).neighbor_front;
                            PHTelem(right_front_neighbors(2)).C = bezt.Ce4;
                            PHTelem(right_front_neighbors(4)).C = bezt.Ce8;
                            PHTelem(right_front_neighbors(2)).nodes(reg.ne_mid) = newBasisSet{4}(reg.ne_mid);
                            PHTelem(right_front_neighbors(2)).nodes(reg.ne_top) = newBasisSet{4}(reg.ne_top);
                            PHTelem(right_front_neighbors(4)).nodes(reg.ne_low) = newBasisSet{8}(reg.ne_low);
                            PHTelem(right_front_neighbors(4)).nodes(reg.ne_mid) = newBasisSet{8}(reg.ne_mid);
                        end
                    case 18
                        if ~isempty(PHTelem(e).neighbor_left_back) && (parentIndex==1)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 17, elmIn_n, p, q, r, dimBasisTemp, reg );
                            left_back_neighbors = PHTelem(e).neighbor_left_back;
                            PHTelem(left_back_neighbors(1)).C = bezt.Ce2;
                            PHTelem(left_back_neighbors(2)).C = bezt.Ce6;
                            PHTelem(left_back_neighbors(1)).nodes(reg.se_mid) = newBasisSet{2}(reg.se_mid);
                            PHTelem(left_back_neighbors(1)).nodes(reg.se_top) = newBasisSet{2}(reg.se_top);
                            PHTelem(left_back_neighbors(2)).nodes(reg.se_low) = newBasisSet{6}(reg.se_low);
                            PHTelem(left_back_neighbors(2)).nodes(reg.se_mid) = newBasisSet{6}(reg.se_mid);
                        end
                        if ~isempty(PHTelem(e).neighbor_left) && (parentIndex==2)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 19, elmIn_n, p, q, r, dimBasisTemp, reg );
                            left_back_neighbors = PHTelem(e).neighbor_left;
                            PHTelem(left_back_neighbors(2)).C = bezt.Ce4;
                            PHTelem(left_back_neighbors(4)).C = bezt.Ce8;
                            PHTelem(left_back_neighbors(2)).nodes(reg.ne_mid) = newBasisSet{4}(reg.ne_mid);
                            PHTelem(left_back_neighbors(2)).nodes(reg.ne_top) = newBasisSet{4}(reg.ne_top);
                            PHTelem(left_back_neighbors(4)).nodes(reg.ne_low) = newBasisSet{8}(reg.ne_low);
                            PHTelem(left_back_neighbors(4)).nodes(reg.ne_mid) = newBasisSet{8}(reg.ne_mid);
                        end
                        if ~isempty(PHTelem(e).neighbor_back) && (parentIndex==3)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 16, elmIn_n, p, q, r, dimBasisTemp, reg );
                            left_back_neighbors = PHTelem(e).neighbor_back;
                            PHTelem(left_back_neighbors(1)).C = bezt.Ce1;
                            PHTelem(left_back_neighbors(3)).C = bezt.Ce5;
                            PHTelem(left_back_neighbors(1)).nodes(reg.sw_mid) = newBasisSet{1}(reg.sw_mid);
                            PHTelem(left_back_neighbors(1)).nodes(reg.sw_top) = newBasisSet{1}(reg.sw_top);
                            PHTelem(left_back_neighbors(3)).nodes(reg.sw_low) = newBasisSet{5}(reg.sw_low);
                            PHTelem(left_back_neighbors(3)).nodes(reg.sw_mid) = newBasisSet{5}(reg.sw_mid);
                        end
                    case 19
                        if ~isempty(PHTelem(e).neighbor_right_back) && (parentIndex==1)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 16, elmIn_n, p, q, r, dimBasisTemp, reg );
                            right_back_neighbors = PHTelem(e).neighbor_right_back;
                            PHTelem(right_back_neighbors(1)).C = bezt.Ce1;
                            PHTelem(right_back_neighbors(2)).C = bezt.Ce5;
                            PHTelem(right_back_neighbors(1)).nodes(reg.sw_mid) = newBasisSet{1}(reg.sw_mid);
                            PHTelem(right_back_neighbors(1)).nodes(reg.sw_top) = newBasisSet{1}(reg.sw_top);
                            PHTelem(right_back_neighbors(2)).nodes(reg.sw_low) = newBasisSet{5}(reg.sw_low);
                            PHTelem(right_back_neighbors(2)).nodes(reg.sw_mid) = newBasisSet{5}(reg.sw_mid);
                        end
                        if ~isempty(PHTelem(e).neighbor_right) && (parentIndex==2)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 18, elmIn_n, p, q, r, dimBasisTemp, reg );
                            right_back_neighbors = PHTelem(e).neighbor_right;
                            PHTelem(right_back_neighbors(2)).C = bezt.Ce3;
                            PHTelem(right_back_neighbors(4)).C = bezt.Ce7;
                            PHTelem(right_back_neighbors(2)).nodes(reg.nw_mid) = newBasisSet{3}(reg.nw_mid);
                            PHTelem(right_back_neighbors(2)).nodes(reg.nw_top) = newBasisSet{3}(reg.nw_top);
                            PHTelem(right_back_neighbors(4)).nodes(reg.nw_low) = newBasisSet{7}(reg.nw_low);
                            PHTelem(right_back_neighbors(4)).nodes(reg.nw_mid) = newBasisSet{7}(reg.nw_mid);
                        end
                        if ~isempty(PHTelem(e).neighbor_back) && (parentIndex==3)
                            [newBasisSet, dimBasisTemp2] = newBasisIndices3D( 17, elmIn_n, p, q, r, dimBasisTemp, reg );
                            right_back_neighbors = PHTelem(e).neighbor_back;
                            PHTelem(right_back_neighbors(2)).C = bezt.Ce2;
                            PHTelem(right_back_neighbors(4)).C = bezt.Ce6;
                            PHTelem(right_back_neighbors(2)).nodes(reg.se_mid) = newBasisSet{2}(reg.se_mid);
                            PHTelem(right_back_neighbors(2)).nodes(reg.se_top) = newBasisSet{2}(reg.se_top);
                            PHTelem(right_back_neighbors(4)).nodes(reg.se_low) = newBasisSet{6}(reg.se_low);
                            PHTelem(right_back_neighbors(4)).nodes(reg.se_mid) = newBasisSet{6}(reg.se_mid);
                        end
                end
            end
        end
        dimBasisTemp = dimBasisTemp2;
    end
    
    %calculate the new Bezier extraction operators and element_nod indices
    %of the children elements
    C_temp = PHTelem(e).C;
    elmIn = PHTelem(e).nodes;
    
    [bez,nodes,dimBasis,newBasisIndex] = deCasteljau3dai_coordinates(C_temp,xmin,xmax,ymin,ymax,zmin,zmax,p,q,r,[RTjunct,newBasisVert,7],knotUl,knotUr,knotVd,knotVu,knotWb,knotWt,elmIn,dimBasis);
    indexCounter = 0;
    newBasisIndex = newBasisIndex(:,[1, 3, 2, 4, 5, 7, 6, 8]);
    for indexNewVert = [newBasisVert,RTjunct,7]
        indexCounter = indexCounter + 1;
        coordinates(newBasisIndex(indexNewVert,:),:) = cpts{indexCounter};        
    end
    
    PHTelem(lastElem+1).C = bez.Ce1;
    PHTelem(lastElem+2).C = bez.Ce2;
    PHTelem(lastElem+3).C = bez.Ce3;
    PHTelem(lastElem+4).C = bez.Ce4;
    PHTelem(lastElem+5).C = bez.Ce5;
    PHTelem(lastElem+6).C = bez.Ce6;
    PHTelem(lastElem+7).C = bez.Ce7;
    PHTelem(lastElem+8).C = bez.Ce8;
    PHTelem(lastElem+1).nodes = nodes.in1;
    PHTelem(lastElem+2).nodes = nodes.in2;
    PHTelem(lastElem+3).nodes = nodes.in3;
    PHTelem(lastElem+4).nodes = nodes.in4;
    PHTelem(lastElem+5).nodes = nodes.in5;
    PHTelem(lastElem+6).nodes = nodes.in6;
    PHTelem(lastElem+7).nodes = nodes.in7;
    PHTelem(lastElem+8).nodes = nodes.in8;
    
    %update the neighbors of the neighbors with self
    for ichild = 1:8
        if (length(PHTelem(lastElem+ichild).neighbor_down)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_down).neighbor_up = setdiff([PHTelem(PHTelem(lastElem+ichild).neighbor_down).neighbor_up, lastElem+ichild],e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_right)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_right).neighbor_left = setdiff([PHTelem(PHTelem(lastElem+ichild).neighbor_right).neighbor_left, lastElem+ichild],e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_up)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_up).neighbor_down = setdiff([PHTelem(PHTelem(lastElem+ichild).neighbor_up).neighbor_down, lastElem+ichild],e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_left)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_left).neighbor_right = setdiff([PHTelem(PHTelem(lastElem+ichild).neighbor_left).neighbor_right, lastElem+ichild],e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_front)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_front).neighbor_back = setdiff([PHTelem(PHTelem(lastElem+ichild).neighbor_front).neighbor_back, lastElem+ichild],e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_back)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_back).neighbor_front = setdiff([PHTelem(PHTelem(lastElem+ichild).neighbor_back).neighbor_front, lastElem+ichild],e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_up_left)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_up_left).neighbor_down_right = setdiff([PHTelem(PHTelem(lastElem+ichild).neighbor_up_left).neighbor_down_right, lastElem+ichild],e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_down_left)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_down_left).neighbor_up_right = setdiff([PHTelem(PHTelem(lastElem+ichild).neighbor_down_left).neighbor_up_right, lastElem+ichild],e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_up_right)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_up_right).neighbor_down_left = setdiff([PHTelem(PHTelem(lastElem+ichild).neighbor_up_right).neighbor_down_left, lastElem+ichild],e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_down_right)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_down_right).neighbor_up_left = setdiff([PHTelem(PHTelem(lastElem+ichild).neighbor_down_right).neighbor_up_left, lastElem+ichild],e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_up_front)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_up_front).neighbor_down_back = setdiff([PHTelem(PHTelem(lastElem+ichild).neighbor_up_front).neighbor_down_back, lastElem+ichild],e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_down_front)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_down_front).neighbor_up_back = setdiff([PHTelem(PHTelem(lastElem+ichild).neighbor_down_front).neighbor_up_back, lastElem+ichild],e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_up_back)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_up_back).neighbor_down_front = setdiff([PHTelem(PHTelem(lastElem+ichild).neighbor_up_back).neighbor_down_front, lastElem+ichild],e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_down_back)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_down_back).neighbor_up_front = setdiff([PHTelem(PHTelem(lastElem+ichild).neighbor_down_back).neighbor_up_front, lastElem+ichild],e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_left_front)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_left_front).neighbor_right_back = setdiff([PHTelem(PHTelem(lastElem+ichild).neighbor_left_front).neighbor_right_back, lastElem+ichild],e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_right_front)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_right_front).neighbor_left_back = setdiff([PHTelem(PHTelem(lastElem+ichild).neighbor_right_front).neighbor_left_back, lastElem+ichild],e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_left_back)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_left_back).neighbor_right_front = setdiff(unique([PHTelem(PHTelem(lastElem+ichild).neighbor_left_back).neighbor_right_front, lastElem+ichild]),e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_right_back)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_right_back).neighbor_left_front = setdiff(unique([PHTelem(PHTelem(lastElem+ichild).neighbor_right_back).neighbor_left_front, lastElem+ichild]),e);
        end
    end
end
%TODO: Change if using rational functions
controlPts = coordinates(:,1:4);
fieldData = coordinates(:,5:end);