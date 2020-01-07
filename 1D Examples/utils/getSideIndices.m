function [west,center,east] = getSideIndices(p)
% Calculate regions of a (p+1)*1 matrix

indexMatrix = 1:p+1;
alpha = floor((p-1)/2);
% Pick the sides
west = indexMatrix(1:alpha+1);
center = indexMatrix((alpha+2):p-alpha);
east = indexMatrix(p-alpha+1:p+1);

end
