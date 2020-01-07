function [dirichlet] = initialBC_tensile(PHTelem,geometry)
% Initialize the crack and boundary conditions

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
            bottomEdge = [bottomEdge, PHTelem{1}(i).nodes(down_nodes)];
        end
    end
end

% Set the top boundary degrees of freedom
for i=1:length(PHTelem{1})
    if isempty(PHTelem{1}(i).children)
        if isempty(PHTelem{1}(i).neighbor_up)
            topEdge = [topEdge, PHTelem{1}(i).nodes(up_nodes)];
        end
    end
end

% Set the left boundary degrees of freedom
for i=1:length(PHTelem{1})
    if isempty(PHTelem{1}(i).children)
        if isempty(PHTelem{1}(i).neighbor_left)
            leftEdge = [leftEdge, PHTelem{1}(i).nodes(left_nodes)];
        end
    end
end

% Set the right boundary degrees of freedom
for i=1:length(PHTelem{1})
    if isempty(PHTelem{1}(i).children)
        if isempty(PHTelem{1}(i).neighbor_right)
            rightEdge = [rightEdge, PHTelem{1}(i).nodes(right_nodes)];
        end
    end
end
dirichlet.Left = unique(leftEdge,'stable');
dirichlet.Right = unique(rightEdge,'stable');
dirichlet.Top = unique(topEdge,'stable');
dirichlet.Bottom = unique(bottomEdge,'stable');

dirichlet.XY =[dirichlet.Bottom,dirichlet.Top,dirichlet.Left]; % Contains the node numbers whose DOF/s have to be restrained.
dirichlet.ValXY = zeros(length(dirichlet.XY'),2); %  0 = Unrestrained 1 = Restrained. X and Y DOFs of the nodes are considered.

dirichlet.ValXY (1:2*length(dirichlet.Bottom),2) = 1; % Restrained Y
dirichlet.ValXY (2*length(dirichlet.Bottom)+1:end,1) = 1; % Restrained X

dirichlet.restrainedPts = zeros(size(dirichlet.ValXY));% Prescribed values of Contrained Dofs.
dirichlet.restrainedPts(length(dirichlet.Bottom)+1:2*length(dirichlet.Bottom),2) =1;
dirichlet.reactForce = dirichlet.restrainedPts;

end