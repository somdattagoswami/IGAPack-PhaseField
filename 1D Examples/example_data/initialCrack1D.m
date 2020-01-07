function [dirichlet] = initialCrack1D(PHTelem)
%initialize the crack and boundary conditions

centre = [];
% Set the bottom boundary degrees of freedom
for i=1:length(PHTelem)
    if isempty(PHTelem(i).children)
        if isempty(PHTelem(i).neighbor_left)
            left = PHTelem(i).nodes(1);
            centre = [centre,PHTelem(i).nodes(2:end);];
        else
            centre = [centre,PHTelem(i).nodes;];
        end
        
    end
end

% Set the left boundary degrees of freedom
for i=1:length(PHTelem)
    if isempty(PHTelem(i).children)
        if isempty(PHTelem(i).neighbor_right)
            right = PHTelem(i).nodes(end);
        end
    end
end

dirichlet.Left = left;
dirichlet.Right = right;
dirichlet.Centre = unique(centre(1:end-1),'stable');

dirichlet.X =[dirichlet.Left,dirichlet.Right]; % Contains the node numbers whose DOF/s have to be restrained.
dirichlet.ValX = zeros(size(dirichlet.X)); 
end