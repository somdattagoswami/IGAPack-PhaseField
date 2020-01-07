function [dirichlet] = initialBC_shear(PHTelem,geometry)
% Initialize the boundary conditions
p = geometry.p;
q = geometry.q;
% Define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

bottomEdge = [];
topEdge = [];
leftEdge = [];
rightEdge = [];

% Set the bottom boundary degrees of freedom

for i=1:length(PHTelem{1})
    if isempty(PHTelem{1}(i).children)
        if isempty(PHTelem{1}(i).neighbor_down)
            bottomEdge = [bottomEdge, PHTelem{1}(i).nodesGlobal(down_nodes)];
        end
        if isempty(PHTelem{1}(i).neighbor_up)
            topEdge = [topEdge, PHTelem{1}(i).nodesGlobal(up_nodes)];
        end
        if isempty(PHTelem{1}(i).neighbor_left)
            leftEdge = [leftEdge, PHTelem{1}(i).nodesGlobal(left_nodes)];
        end
        if isempty(PHTelem{1}(i).neighbor_right)
            rightEdge = [rightEdge, PHTelem{1}(i).nodesGlobal(right_nodes)];
        end
    end
end



dirichlet.Left = unique(leftEdge,'stable');
dirichlet.Right = unique(rightEdge,'stable');
dirichlet.Top = unique(topEdge,'stable');
dirichlet.Bottom = unique(bottomEdge,'stable');


dirichlet.XY =[dirichlet.Bottom,dirichlet.Top,dirichlet.Left,dirichlet.Right]; % Contains the node numbers whose DOF/s have to be restrained.
dirichlet.ValXY = zeros(length(dirichlet.XY'),2); %  0 = Unrestrained 1 = Restrained. X and Y DOFs of the nodes are considered.

dirichlet.ValXY (1:length(dirichlet.Bottom)+length(dirichlet.Top),1) = 1; % Restrained X
dirichlet.ValXY (1:end,2) = 1; % Restrained Y

dirichlet.restrainedPts = zeros(size(dirichlet.ValXY));% Prescribed values of Contrained Dofs.
dirichlet.restrainedPts(length(dirichlet.Bottom)+1:length(dirichlet.Bottom)+length(dirichlet.Top),1) =1;
dirichlet.reactForce = zeros(size(dirichlet.ValXY));
dirichlet.reactForce(1:length(dirichlet.Bottom),1) = 1;
end