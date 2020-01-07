function [dirichlet] = initialBC_cube(PHTelem,geometry)
% Initialize the crack and boundary conditions

p = geometry.p;
q = geometry.q;
r = geometry.r;

% Defining side node indices
nodesBottom = 1:(p+1)*(q+1);
nodesRight = (p+1):(p+1):(p+1)*(q+1)*(r+1);
nodesTop = 1+(p+1)*(q+1)*r:(p+1)*(q+1)*(r+1);
nodesLeft = 1:(p+1):(1+(p+1)*(q+1)*(r+1)-p);
nodesFront = [];

for i=1:r+1
    nodesFront = [nodesFront, (p+1)*(q+1)*(i-1)+1:(p+1)*((q+1)*(i-1)+1)];
end
nodesBack = nodesFront + (p+1)*q;

faceBottom = [];
faceTop = [];
faceLeft = [];
faceRight = [];
faceBack = [];

for patchIndex = 1
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_down)
                faceBottom = [faceBottom, PHTelem{patchIndex}(i).nodesGlobal(nodesBottom)];
            end
            if isempty(PHTelem{patchIndex}(i).neighbor_up)
                faceTop = [faceTop, PHTelem{patchIndex}(i).nodesGlobal(nodesTop)];
            end
            if isempty(PHTelem{patchIndex}(i).neighbor_left)
                faceLeft = [faceLeft, PHTelem{patchIndex}(i).nodesGlobal(nodesLeft)];
            end
             if isempty(PHTelem{patchIndex}(i).neighbor_right)
                faceRight = [faceRight, PHTelem{patchIndex}(i).nodesGlobal(nodesRight)];
             end
            if isempty(PHTelem{patchIndex}(i).neighbor_back)
                faceBack = [faceBack, PHTelem{patchIndex}(i).nodesGlobal(nodesBack)];
            end
        end
    end
end

dirichlet.Left = unique(faceLeft,'stable');
dirichlet.Right = unique(faceRight,'stable');
dirichlet.Top = unique(faceTop,'stable');
dirichlet.Bottom = unique(faceBottom,'stable');
dirichlet.Back = unique(faceBack,'stable');

dirichlet.XYZ =[dirichlet.Bottom,dirichlet.Top,dirichlet.Back,dirichlet.Left]; % Contains the node numbers whose DOF/s have to be restrained.
dirichlet.ValXYZ = zeros(length(dirichlet.XYZ),3); %  0 = Unrestrained 1 = Restrained. X and Y DOFs of the nodes are considered.
dirichlet.ValXYZ(1:2*length(dirichlet.Bottom),3) = 1; % Restrained Z
dirichlet.ValXYZ(2*length(dirichlet.Bottom)+1:2*length(dirichlet.Bottom)+length(dirichlet.Back),2) = 1; % Restrained Y
dirichlet.ValXYZ (2*length(dirichlet.Bottom)+length(dirichlet.Back)+1:end,1) = 1; % Restrained X
dirichlet.restrainedPts = zeros(size(dirichlet.ValXYZ));% Prescribed values of Contrained Dofs.
dirichlet.restrainedPts(length(dirichlet.Bottom)+1:2*length(dirichlet.Bottom),3) =1;
dirichlet.reactForce = dirichlet.restrainedPts;
end