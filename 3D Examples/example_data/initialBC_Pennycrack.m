function [dirichlet] = initialBC_Pennycrack(PHTelem,geometry)
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

for patchIndex = 1
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_down)
                faceBottom = [faceBottom, PHTelem{patchIndex}(i).nodesGlobal(nodesBottom)];
            end
            if isempty(PHTelem{patchIndex}(i).neighbor_up)
                faceTop = [faceTop, PHTelem{patchIndex}(i).nodesGlobal(nodesTop)];
            end
        end
    end
end

dirichlet.Top = unique(faceTop,'stable');
dirichlet.Bottom = unique(faceBottom,'stable');

dirichlet.XYZ =[dirichlet.Bottom,dirichlet.Top]; % Contains the node numbers whose DOF/s have to be restrained.
dirichlet.ValXYZ = ones(length(dirichlet.XYZ),3); %  0 = Unrestrained 1 = Restrained.
dirichlet.restrainedPts = zeros(size(dirichlet.ValXYZ));% Prescribed values of Contrained Dofs.
dirichlet.restrainedPts(1:length(dirichlet.Bottom),3) = -1;
dirichlet.restrainedPts(length(dirichlet.Bottom)+1:end,3) = 1;
dirichlet.reactForce = dirichlet.restrainedPts;
end