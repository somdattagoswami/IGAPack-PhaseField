function [sw,south,se,west,center,east,nw,north,ne] = getCornerIndices(p,q)

%calculate regions of a (p+1)*(q+1) matrix

alpha = floor((p-1)/2); 
beta = floor((q-1)/2);


%form the index matrix
%indexMatrix = reshape(1:(p+alpha)*(q+beta),q+alpha,p+beta)';

indexMatrix = reshape(1:(p+1)*(q+1),q+1,p+1)';

%pick the corners
sw = indexMatrix(1:alpha+1,1:alpha+1);
nw = indexMatrix(q-alpha+1:q+1,1:alpha+1);
ne = indexMatrix(q-alpha+1:q+1,p-alpha+1:p+1);
se = indexMatrix(1:alpha+1,p-alpha+1:p+1);
south = indexMatrix(1:alpha+1, (alpha+2):(p-alpha));
west = indexMatrix(alpha+2:q-alpha, 1:alpha+1);
center = indexMatrix(alpha+2:q-alpha, (alpha+2):p-alpha);
east = indexMatrix(alpha+2:q-alpha,p-alpha+1:p+1);
north = indexMatrix(q-alpha+1:q+1, alpha+2:p-alpha);
% 
% 
% %convert to row vectors
sw = reshape(sw,1,(alpha+1)*(beta+1));
se = reshape(se,1,(alpha+1)*(beta+1));
ne = reshape(ne,1,(alpha+1)*(beta+1));
nw = reshape(nw,1,(alpha+1)*(beta+1));
north = reshape(north,1,(beta+1)*(p-2*alpha-1));
west = reshape(west, 1, (beta+1)*(p-2*alpha-1));
south = reshape(south, 1, (beta+1)*(p-2*alpha-1));
east = reshape(east, 1, (beta+1)*(p-2*alpha-1));


end

