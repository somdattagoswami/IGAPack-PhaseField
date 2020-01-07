function [dirichlet] = initialBC_cantilever(PHTelem,geometry,controlPts)
% Initialize the crack and boundary conditions

p = geometry.p;
q = geometry.q;
% Define side node indices
right_nodes = (p+1):(p+1):(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

leftTop = [];
leftBottom = [];
support = [];
% Set the bottom boundary degrees of freedom
for i=1:length(PHTelem{1})
    if isempty(PHTelem{1}(i).children)
        if isempty(PHTelem{1}(i).neighbor_right)
            rightEdge = PHTelem{1}(i).nodesGlobal(right_nodes);
            for iPoint = 1:length(rightEdge)
                if controlPts{1}(rightEdge(iPoint),2) > 0.40 && controlPts{1}(rightEdge(iPoint),2)< 0.6
                support = [support,rightEdge(iPoint)];
                end
            end
        end
    end
end

for i=1:length(PHTelem{1})
    if isempty(PHTelem{1}(i).children)
        if isempty(PHTelem{1}(i).neighbor_left)
            leftEdge = PHTelem{1}(i).nodesGlobal(left_nodes);
            for iPoint = 1:length(leftEdge)
                controlPts{1}(leftEdge(iPoint),2);
                if controlPts{1}(leftEdge(iPoint),2) >= 0.5
                    leftTop = [leftTop, leftEdge(iPoint)];
                else
                    leftBottom = [leftBottom, leftEdge(iPoint)];
                end
            end
        end
    end
end
dirichlet.leftTop = unique(leftTop,'stable');
dirichlet.leftBottom = unique(leftBottom,'stable');

dirichlet.XY =[support,dirichlet.leftTop,dirichlet.leftBottom]; % Contains the node numbers whose DOF/s have to be restrained.
dirichlet.ValXY = zeros(length(dirichlet.XY'),2); %  0 = Unrestrained 1 = Restrained. X and Y DOFs of the nodes are considered.
dirichlet.ValXY (1,1:2) = 1; % Restrained X and Y
dirichlet.ValXY (length(support)+1:end,2) = 1; % Restrained Y

dirichlet.restrainedPts = zeros(size(dirichlet.ValXY));% Prescribed values of Contrained Dofs.
dirichlet.restrainedPts(length(support)+1:length(support)+length(dirichlet.leftTop),2) = 1;
dirichlet.restrainedPts(length(support)+length(dirichlet.leftTop)+1:end,2) = -1;
dirichlet.reactForce = zeros(size(dirichlet.ValXY));
dirichlet.reactForce = dirichlet.restrainedPts;
end