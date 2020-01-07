 function [PHTelem,controlPts,dimBasis,octupleList] = genControlPts3D(nurbs,geometry)
% Generates the controlPts array and the knot vectors for the C1 mesh
% No initial repeating of control points
% No degree elevation

dim = 3; % Number of physical dimensions
p = geometry.p;
q = geometry.q;
r = geometry.r;
toleq = 1e-10;% Tolerance for equality tests

knotU = nurbs.knots{1};
knotV = nurbs.knots{2};
knotW = nurbs.knots{3};
coefs = nurbs.coefs;

% The number of control points in the u,v and w directions
basisU = length(knotU)-p-1; % Number of basis functions in the u direction
basisV = length(knotV)-q-1; % Number of basis functions in the v direction
basisW = length(knotW)-r-1; % Number of basis functions in the w direction
dimBasis = basisU*basisV*basisW; % Total number of control points
controlPts = zeros(dimBasis,dim+1); %allocate one more column for the weights
numElements = geometry.numElmtU*geometry.numElmtV*geometry.numElmtW;
octupleList = zeros(numElements,8);

index = 0;
for k=1:basisW % For each node in the z direction
    for j=1:basisV % For each node in the y direction
        for i=1:basisU % For each node in the x direction
            index = index + 1; % The index of the coordinate array
            controlPts(index,:) = [coefs(1,i,j,k)./coefs(4,i,j,k), coefs(2,i,j,k)./coefs(4,i,j,k), coefs(3,i,j,k)./coefs(4,i,j,k), coefs(4,i,j,k)]; 
        end % For i
    end % For j
end % For k

% Do Bezier extraction
[C_u, ~] = bezierExtraction(knotU,geometry.p);
[C_v, ~] = bezierExtraction(knotV,geometry.q);
[C_w, ~] = bezierExtraction(knotW,geometry.r);

PHTelem = struct;
% Initialize the neighbor connectivity lists
PHTelem.neighbor_left = [];
PHTelem.neighbor_right = [];
PHTelem.neighbor_down = [];
PHTelem.neighbor_up = [];
PHTelem.neighbor_front = [];
PHTelem.neighbor_back = [];
PHTelem.neighbor_up_left = [];
PHTelem.neighbor_down_left = [];
PHTelem.neighbor_up_right = [];
PHTelem.neighbor_down_right = [];
PHTelem.neighbor_up_front = [];
PHTelem.neighbor_down_front = [];
PHTelem.neighbor_up_back = [];
PHTelem.neighbor_down_back = [];
PHTelem.neighbor_left_front = [];
PHTelem.neighbor_right_front = [];
PHTelem.neighbor_left_back = [];
PHTelem.neighbor_right_back = [];


% Loop through each element and compute the element-node connectivities
nument = (p+1)*(q+1)*(r+1);
elementCounter = 0;
for k=1:length(knotW)-1
    for j=1:length(knotV)-1
        for i=1:length(knotU)-1
            if (knotU(i+1)>knotU(i)) && (knotV(j+1) >knotV(j)) && (knotW(k+1) > knotW(k))  %the knotspan has non-zero area
                elementCounter = elementCounter + 1;
                PHTelem(elementCounter).parent = [];
                PHTelem(elementCounter).children = [];
                PHTelem(elementCounter).vertex = [knotU(i), knotV(j), knotW(k), knotU(i+1), knotV(j+1), knotW(k+1)];
                
                tcount = 0;
                currow = zeros(1, nument);
                % Now we add the nodes from i-p...i in the u
                % direction, j-q...j in the v direction, k-r...k in
                % w direction
                for t3 = k-r:k
                    for t2=j-q:j
                        for t1 = i-p:i
                            tcount = tcount + 1;
                            currow(tcount) = t1+(t2-1)*basisU+(t3-1)*basisU*basisV;
                        end
                    end
                end
                PHTelem(elementCounter).nodes=currow;
                PHTelem(elementCounter).level = 0;
                octupleList(elementCounter,:) = numElements + (8*(elementCounter-1)+1:8*elementCounter);
            end
        end
    end
end


% Loop through each element and compute the neighbor lists and Bezier extraction operators
indexMatrix = permute(reshape(1:numElements,geometry.numElmtU,geometry.numElmtV,geometry.numElmtW),[2,1,3]);
for k=1:geometry.numElmtW
    for j=1:geometry.numElmtV
        for i=1:geometry.numElmtU
            elementIndex = indexMatrix(j,i,k);
            PHTelem(elementIndex).C = kron(kron(C_w(:,:,k),C_v(:,:,j)),C_u(:,:,i));
            
            if i>1
                PHTelem(elementIndex).neighbor_left = indexMatrix(j,i-1,k);
                if j>1
                    PHTelem(elementIndex).neighbor_left_front = indexMatrix(j-1,i-1,k);
                end
                if j<geometry.numElmtV
                    PHTelem(elementIndex).neighbor_left_back = indexMatrix(j+1,i-1,k);
                end
                if k>1
                    PHTelem(elementIndex).neighbor_down_left = indexMatrix(j,i-1,k-1);
                end
                if k<geometry.numElmtW
                    PHTelem(elementIndex).neighbor_up_left = indexMatrix(j,i-1,k+1);
                end
            end
            if i<geometry.numElmtU
                PHTelem(elementIndex).neighbor_right = indexMatrix(j,i+1,k);
                if j>1
                    PHTelem(elementIndex).neighbor_right_front = indexMatrix(j-1,i+1,k);
                end
                if j<geometry.numElmtV
                    PHTelem(elementIndex).neighbor_right_back = indexMatrix(j+1,i+1,k);
                end
                if k>1
                    PHTelem(elementIndex).neighbor_down_right = indexMatrix(j,i+1,k-1);
                end
                if k<geometry.numElmtW
                    PHTelem(elementIndex).neighbor_up_right = indexMatrix(j,i+1,k+1);
                end
            end
            
            if k>1
                PHTelem(elementIndex).neighbor_down = indexMatrix(j,i,k-1);
                if j>1
                    PHTelem(elementIndex).neighbor_down_front = indexMatrix(j-1,i,k-1);
                end
                if j<geometry.numElmtV
                    PHTelem(elementIndex).neighbor_down_back = indexMatrix(j+1,i,k-1);
                end
            end
            
            if k<geometry.numElmtW
                PHTelem(elementIndex).neighbor_up = indexMatrix(j,i,k+1);
                if j>1
                    PHTelem(elementIndex).neighbor_up_front = indexMatrix(j-1,i,k+1);
                end
                if j<geometry.numElmtV
                    PHTelem(elementIndex).neighbor_up_back = indexMatrix(j+1,i,k+1);
                end
            end
            
            if j>1
                PHTelem(elementIndex).neighbor_front = indexMatrix(j-1,i,k);
            end
            if j<geometry.numElmtV
                PHTelem(elementIndex).neighbor_back = indexMatrix(j+1,i,k);
            end
        end
    end
end

[PHTelem,dimBasis,controlPts] = crossInsert3D(PHTelem,controlPts,1:numElements,dimBasis,geometry);

end


