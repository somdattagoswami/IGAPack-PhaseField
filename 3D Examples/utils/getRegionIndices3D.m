function [reg] = getRegionIndices3D(p,q,r)
%calculate regions of a (p+1)*(q+1)*(r+1) matrix

reg = struct;

alpha = floor((p-1)/2);
beta = floor((q-1)/2);
gamma = floor((r-1)/2);

%form the index matrix
indexMatrix = permute(reshape(1:(p+1)*(q+1)*(r+1),r+1,q+1,p+1),[2,1,3]);

%pick the corners on the lower slice
reg.sw_low = indexMatrix(1:alpha+1,1:beta+1,1:gamma+1);
reg.nw_low = indexMatrix(q-beta+1:q+1,1:beta+1,1:alpha+1);
reg.ne_low = indexMatrix(q-beta+1:q+1,p-alpha+1:p+1,1:alpha+1);
reg.se_low = indexMatrix(1:alpha+1,p-alpha+1:p+1,1:alpha+1);

%pick the center regions on the lower slice
reg.south_low = indexMatrix(1:alpha+1, alpha+2:p-alpha,1:alpha+1);
reg.west_low = indexMatrix(alpha+2:q-beta, 1:alpha+1,1:alpha+1);
reg.center_low = indexMatrix(alpha+2:q-beta, alpha+2:p-alpha,1:alpha+1);
reg.east_low = indexMatrix(alpha+2:q-beta,p-alpha+1:p+1,1:alpha+1);
reg.north_low = indexMatrix(q-beta+1:q+1, alpha+2:p-alpha,1:alpha+1);

%pick the corners on the middle slice
reg.sw_mid = indexMatrix(1:alpha+1,1:alpha+1,gamma+2:r-gamma);
reg.nw_mid = indexMatrix(q-beta+1:q+1,1:alpha+1,gamma+2:r-gamma);
reg.ne_mid = indexMatrix(q-beta+1:q+1,p-alpha+1:p+1,gamma+2:r-gamma);
reg.se_mid = indexMatrix(1:alpha+1,p-alpha+1:p+1,gamma+2:r-gamma);

%pick the center regions on the middle slice
reg.south_mid = indexMatrix(1:alpha+1, beta+2:p-beta,gamma+2:r-gamma);
reg.west_mid = indexMatrix(beta+2:q-beta, 1:alpha+1,gamma+2:r-gamma);
reg.center_mid = indexMatrix(beta+2:q-beta, alpha+2:p-alpha,gamma+2:r-gamma);
reg.east_mid = indexMatrix(beta+2:q-beta,p-alpha+1:p+1,gamma+2:r-gamma);
reg.north_mid = indexMatrix(q-beta+1:q+1, alpha+2:p-alpha,gamma+2:r-gamma);

%pick the corners on the top slice
reg.sw_top = indexMatrix(1:alpha+1,1:beta+1,r-gamma+1:r+1);
reg.nw_top = indexMatrix(q-beta+1:q+1,1:beta+1,r-gamma+1:r+1);
reg.ne_top = indexMatrix(q-beta+1:q+1,p-alpha+1:p+1,r-gamma+1:r+1);
reg.se_top = indexMatrix(1:alpha+1,p-beta+1:p+1,r-gamma+1:r+1);

%pick the center regions on the top slice
reg.south_top = indexMatrix(1:alpha+1, alpha+2:p-alpha,r-gamma+1:r+1);
reg.west_top = indexMatrix(beta+2:q-beta, 1:beta+1,r-gamma+1:r+1);
reg.center_top = indexMatrix(beta+2:q-beta, alpha+2:p-alpha,r-gamma+1:r+1);
reg.east_top = indexMatrix(beta+2:q-beta,p-alpha+1:p+1,r-gamma+1:r+1);
reg.north_top = indexMatrix(q-beta+1:q+1, alpha+2:p-alpha,r-gamma+1:r+1);

%convert to row vectors
reg.sw_low = reshape(reg.sw_low,1,(alpha+1)*(beta+1)*(gamma+1));
reg.se_low = reshape(reg.se_low,1,(alpha+1)*(beta+1)*(gamma+1));
reg.ne_low = reshape(reg.ne_low,1,(alpha+1)*(beta+1)*(gamma+1));
reg.nw_low = reshape(reg.nw_low,1,(alpha+1)*(beta+1)*(gamma+1));
reg.north_low = reshape(reg.north_low,1,(beta+1)*(gamma+1)*(p-2*alpha-1));
reg.west_low = reshape(reg.west_low, 1, (alpha+1)*(gamma+1)*(q-2*beta-1));
reg.south_low = reshape(reg.south_low, 1, (alpha+1)*(gamma+1)*(p-2*alpha-1));
reg.east_low = reshape(reg.east_low, 1, (alpha+1)*(gamma+1)*(q-2*beta-1));
reg.center_low = reshape(reg.center_low, 1, (gamma+1)*(p-2*alpha-1)*(q-2*beta-1));

reg.sw_mid = reshape(reg.sw_mid,1,(alpha+1)*(beta+1)*(r-2*gamma-1));
reg.se_mid = reshape(reg.se_mid,1,(alpha+1)*(beta+1)*(r-2*gamma-1));
reg.ne_mid = reshape(reg.ne_mid,1,(alpha+1)*(beta+1)*(r-2*gamma-1));
reg.nw_mid = reshape(reg.nw_mid,1,(alpha+1)*(beta+1)*(r-2*gamma-1));
reg.north_mid = reshape(reg.north_mid,1,(beta+1)*(p-2*alpha-1)*(r-2*gamma-1));
reg.west_mid = reshape(reg.west_mid, 1, (alpha+1)*(q-2*beta-1)*(r-2*gamma-1));
reg.south_mid = reshape(reg.south_mid, 1, (beta+1)*(p-2*alpha-1)*(r-2*gamma-1));
reg.east_mid = reshape(reg.east_mid, 1, (alpha+1)*(q-2*beta-1)*(r-2*gamma-1));
reg.center_mid = reshape(reg.center_mid, 1, (p-2*alpha-1)*(q-2*beta-1)*(r-2*gamma-1));

reg.sw_top = reshape(reg.sw_top,1,(alpha+1)*(beta+1)*(gamma+1));
reg.se_top = reshape(reg.se_top,1,(alpha+1)*(beta+1)*(gamma+1));
reg.ne_top = reshape(reg.ne_top,1,(alpha+1)*(beta+1)*(gamma+1));
reg.nw_top = reshape(reg.nw_top,1,(alpha+1)*(beta+1)*(gamma+1));
reg.north_top = reshape(reg.north_top,1,(beta+1)*(gamma+1)*(p-2*alpha-1));
reg.west_top = reshape(reg.west_top, 1, (alpha+1)*(gamma+1)*(q-2*beta-1));
reg.south_top = reshape(reg.south_top, 1, (beta+1)*(gamma+1)*(p-2*alpha-1));
reg.east_top = reshape(reg.east_top, 1, (alpha+1)*(gamma+1)*(q-2*beta-1));
reg.center_top = reshape(reg.center_top, 1, (p-2*alpha-1)*(q-2*beta-1)*(gamma+1));


end

