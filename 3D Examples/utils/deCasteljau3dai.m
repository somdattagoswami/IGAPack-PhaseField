function [bez,nodes,dimBasis] = deCasteljau3dai(Ce,knotU1,knotU2,knotV1,knotV2,knotW1,knotW2,p,q,r,newVert,knotUl,knotUr,knotVd,knotVu,knotWb,knotWt,elmIn,dimBasis)
% Split element with Bezier extraction operator Ce into 8 elements with
% Bezier extraction operator bez.Ce1, bez.Ce2, ..., bez.Ce8
%newVert: new Vertices
%automatic detection of extra basis vertices (saved entries)
%calculates element_nod indices for new elements
%encoding scheme newVert: 1-front, 2-right, 3-back, 4-left, 5-down, 6-up, 7-center
%8-up_left, 9-down_left,


num_rows = size(Ce,1);
num_cols = size(Ce,2);

bez.Ce1 = zeros(num_rows, num_cols);
bez.Ce2 = zeros(num_rows, num_cols);
bez.Ce3 = zeros(num_rows, num_cols);
bez.Ce4 = zeros(num_rows, num_cols);
bez.Ce5 = zeros(num_rows, num_cols);
bez.Ce6 = zeros(num_rows, num_cols);
bez.Ce7 = zeros(num_rows, num_cols);
bez.Ce8 = zeros(num_rows, num_cols);

knotUmid = (knotU1+knotU2)/2;
knotVmid = (knotV1+knotV2)/2;
knotWmid = (knotW1+knotW2)/2;

newKnotU = [knotUl*ones(1,p+1), knotU1*ones(1,p-1), knotUmid*ones(1,p-1), knotU2*ones(1,p-1), knotUr*ones(1, p+1)];
newKnotV = [knotVd*ones(1,q+1), knotV1*ones(1,q-1), knotVmid*ones(1,q-1), knotV2*ones(1,q-1), knotVu*ones(1, q+1)];
newKnotW = [knotWb*ones(1,r+1), knotW1*ones(1,r-1), knotWmid*ones(1,r-1), knotW2*ones(1,r-1), knotWt*ones(1, q+1)];

[newC_u, ~] = bezierExtraction(newKnotU, p);
[newC_v, ~] = bezierExtraction(newKnotV, q);
[newC_w, ~] = bezierExtraction(newKnotW, r);

%calculate the tensor product Bezier extraction operators on the eight
%new elements
newC1 = kron(kron(newC_w(:,:,2),newC_v(:,:,2)),newC_u(:,:,2));
newC2 = kron(kron(newC_w(:,:,2),newC_v(:,:,2)),newC_u(:,:,3));
newC3 = kron(kron(newC_w(:,:,2),newC_v(:,:,3)),newC_u(:,:,2));
newC4 = kron(kron(newC_w(:,:,2),newC_v(:,:,3)),newC_u(:,:,3));
newC5 = kron(kron(newC_w(:,:,3),newC_v(:,:,2)),newC_u(:,:,2));
newC6 = kron(kron(newC_w(:,:,3),newC_v(:,:,2)),newC_u(:,:,3));
newC7 = kron(kron(newC_w(:,:,3),newC_v(:,:,3)),newC_u(:,:,2));
newC8 = kron(kron(newC_w(:,:,3),newC_v(:,:,3)),newC_u(:,:,3));

Ce_save = cell(1,8);

%copy the basis function indices to the children elements
nodes.in1 = elmIn;
nodes.in2 = elmIn;
nodes.in3 = elmIn;
nodes.in4 = elmIn;
nodes.in5 = elmIn;
nodes.in6 = elmIn;
nodes.in7 = elmIn;
nodes.in8 = elmIn;

%do tensor product deCasteljau
for i=1:size(Ce,1)
    
    cur_row = Ce(i,:);     
    [ temp ] = deCasteljau3DTP_mex( cur_row, p, q, r );
    
    %zero out entries coresponding to new basis vertices
    for j=1:length(newVert)
        switch newVert(j)
            case 1
                temp(1:q-1,3:2*p,3:2*r) = 0;
            case 2
                temp(3:2*q,end-p+2:end,3:2*r) = 0;
            case 3
                temp(end-q+2:end,3:2*p,3:2*r) = 0;
            case 4
                temp(3:2*q,1:p-1,3:2*r) = 0;
            case 5
                temp(3:2*q,3:2*p,1:r-1) = 0;
            case 6
                temp(3:2*q,3:2*p,end-r+2:end) = 0;
            case 7
                temp(3:2*q,3:2*p,3:2*r) = 0;
            case 8
                temp(3:2*q,1:p-1,end-r+2:end) = 0;
            case 9
                temp(3:2*q,1:p-1,1:r-1) = 0;
            case 10
                temp(3:2*q,end-p+2:end,end-r+2:end) = 0;
            case 11
                temp(3:2*q,end-p+2:end,1:r-1) = 0;
            case 12
                temp(1:q-1,3:2*p,end-r+2:end) = 0;
            case 13
                temp(1:q-1,3:2*p,1:r-1) = 0;
            case 14
                temp(end-q+2:end,3:2*p,end-r+2:end) = 0;
            case 15
                temp(end-q+2:end,3:2*p,1:r-1) = 0;
            case 16
                temp(1:q-1,1:p-1,3:2*r) = 0;
            case 17
                temp(1:q-1,end-p+2:end,3:2*r) = 0;
            case 18
                temp(end-q+2:end,1:p-1,3:2*r) = 0;
            case 19
                temp(end-q+2:end,end-p+2:end,3:2*r) = 0;
        end
    end
    bez.Ce1(i,:) = reshape(permute(temp(1:q+1,1:p+1,1:r+1),[2,1,3]),1,(p+1)*(q+1)*(r+1));
    bez.Ce2(i,:) = reshape(permute(temp(1:q+1,p+2:end,1:r+1),[2,1,3]),1,(p+1)*(q+1)*(r+1));
    bez.Ce3(i,:) = reshape(permute(temp(q+2:end,1:p+1,1:r+1),[2,1,3]),1,(p+1)*(q+1)*(r+1));
    bez.Ce4(i,:) = reshape(permute(temp(q+2:end,p+2:end,1:r+1),[2,1,3]),1,(p+1)*(q+1)*(r+1));
    bez.Ce5(i,:) = reshape(permute(temp(1:q+1,1:p+1,r+2:end),[2,1,3]),1,(p+1)*(q+1)*(r+1));
    bez.Ce6(i,:) = reshape(permute(temp(1:q+1,p+2:end,r+2:end),[2,1,3]),1,(p+1)*(q+1)*(r+1));
    bez.Ce7(i,:) = reshape(permute(temp(q+2:end,1:p+1,r+2:end),[2,1,3]),1,(p+1)*(q+1)*(r+1));
    bez.Ce8(i,:) = reshape(permute(temp(q+2:end,p+2:end,r+2:end),[2,1,3]),1,(p+1)*(q+1)*(r+1));
    
    Ce_save{1}(i,:) = bez.Ce1(i,:);
    Ce_save{2}(i,:) = bez.Ce2(i,:);
    Ce_save{3}(i,:) = bez.Ce3(i,:);
    Ce_save{4}(i,:) = bez.Ce4(i,:);
    Ce_save{5}(i,:) = bez.Ce5(i,:);
    Ce_save{6}(i,:) = bez.Ce6(i,:);
    Ce_save{7}(i,:) = bez.Ce7(i,:);
    Ce_save{8}(i,:) = bez.Ce8(i,:);
    
end

saveCell = cell(1,8);
savedIndex = cell(1,8);

%add in the new basis functions
%define corners indices
reg = getRegionIndices3D( p,q,r );

[newBasisSet, dimBasis ] = newBasisIndices3D( newVert, elmIn, p, q, r, dimBasis, reg );
%dimBasis

for j=1:length(newVert)
    switch newVert(j)
        case 1
            if norm(bez.Ce1(reg.se_mid,:))>0
                saveCell{1}=[saveCell{1}, 12];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.se_mid)];
            end
            if norm(bez.Ce1(reg.south_mid,:))>0
                saveCell{1}=[saveCell{1}, 11];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.south_mid)];
            end
            if norm(bez.Ce2(reg.south_mid,:))>0
                saveCell{2}=[saveCell{2}, 11];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.south_mid)];
            end
            if norm(bez.Ce2(reg.sw_mid,:))>0
                saveCell{2}=[saveCell{2}, 10];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.sw_mid)];
            end
            if norm(bez.Ce1(reg.se_top,:))>0
                saveCell{1}=[saveCell{1}, 21];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.se_top)];
            end
            if norm(bez.Ce1(reg.south_top,:))>0
                saveCell{1}=[saveCell{1}, 20];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.soth_top)];
            end
            if norm(bez.Ce2(reg.south_top,:))>0
                saveCell{2}=[saveCell{2}, 20];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.south_top)];
            end
            if norm(bez.Ce2(reg.sw_top,:))>0
                saveCell{2}=[saveCell{2}, 19];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.sw_top)];
            end
            if norm(bez.Ce5(reg.se_low,:))>0
                saveCell{5}=[saveCell{5}, 3];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.se_low)];
            end
            if norm(bez.Ce5(reg.south_low,:))>0
                saveCell{5}=[saveCell{5}, 2];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.south_low)];
            end
            if norm(bez.Ce6(reg.south_low,:))>0
                saveCell{6}=[saveCell{6}, 2];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.south_low)];
            end
            if norm(bez.Ce6(reg.sw_low,:))>0
                saveCell{6}=[saveCell{6}, 1];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.sw_low)];
            end
            if norm(bez.Ce5(reg.se_mid,:))>0
                saveCell{5}=[saveCell{5}, 12];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.se_mid)];
            end
            if norm(bez.Ce5(reg.south_mid,:))>0
                saveCell{5}=[saveCell{5}, 11];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.south_mid)];
            end
            if norm(bez.Ce6(reg.south_mid,:))>0
                saveCell{6}=[saveCell{6}, 11];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.south_mid)];
            end
            if norm(bez.Ce6(reg.sw_mid,:))>0
                saveCell{6}=[saveCell{6}, 10];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.sw_mid)];
            end
            
            bez.Ce1(reg.se_mid,:) = newC1(reg.se_mid,:);
            bez.Ce1(reg.south_mid,:) = newC1(reg.south_mid,:);
            bez.Ce2(reg.south_mid,:) = newC2(reg.south_mid,:);
            bez.Ce2(reg.sw_mid,:) = newC2(reg.sw_mid,:);
            bez.Ce1(reg.se_top,:) = newC1(reg.se_top,:);
            bez.Ce1(reg.south_top,:) = newC1(reg.south_top,:);
            bez.Ce2(reg.south_top,:) = newC2(reg.south_top,:);
            bez.Ce2(reg.sw_top,:) = newC2(reg.sw_top,:);
            bez.Ce5(reg.se_low,:) = newC5(reg.se_low,:);
            bez.Ce5(reg.south_low,:) = newC5(reg.south_low,:);
            bez.Ce6(reg.south_low,:) = newC6(reg.south_low,:);
            bez.Ce6(reg.sw_low,:) = newC6(reg.sw_low,:);
            bez.Ce5(reg.se_mid,:) = newC5(reg.se_mid,:);
            bez.Ce5(reg.south_mid,:) = newC5(reg.south_mid,:);
            bez.Ce6(reg.south_mid,:) = newC6(reg.south_mid,:);
            bez.Ce6(reg.sw_mid,:) = newC6(reg.sw_mid,:);
        case 2
            if norm(bez.Ce2(reg.ne_mid,:))>0
                saveCell{2}=[saveCell{2}, 18];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.ne_mid)];
            end
            
            if norm(bez.Ce2(reg.east_mid,:))>0
                saveCell{2}=[saveCell{2}, 15];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.east_mid)];
            end
            
            if norm(bez.Ce4(reg.east_mid,:))>0
                saveCell{4}=[saveCell{4}, 15];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.east_mid)];
            end
            
            if norm(bez.Ce4(reg.se_mid,:))>0
                saveCell{4}=[saveCell{4}, 12];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.se_mid)];
            end
            
            bez.Ce2(reg.ne_mid,:) = newC2(reg.ne_mid,:);
            bez.Ce2(reg.east_mid,:) = newC2(reg.east_mid,:);
            bez.Ce4(reg.east_mid,:) = newC4(reg.east_mid,:);
            bez.Ce4(reg.se_mid,:) = newC4(reg.se_mid,:);
            
            if norm(bez.Ce2(reg.ne_top,:))>0
                saveCell{2}=[saveCell{2}, 27];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.ne_top)];
            end
            
            if norm(bez.Ce2(reg.east_top,:))>0
                saveCell{2}=[saveCell{2}, 24];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.east_top)];
            end
            
            if norm(bez.Ce4(reg.east_top,:))>0
                saveCell{4}=[saveCell{4}, 24];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.east_top)];
            end
            
            if norm(bez.Ce4(reg.se_top,:))>0
                saveCell{4}=[saveCell{4}, 21];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.se_top)];
            end
            
            bez.Ce2(reg.ne_top,:) = newC2(reg.ne_top,:);
            bez.Ce2(reg.east_top,:) = newC2(reg.east_top,:);
            bez.Ce4(reg.east_top,:) = newC4(reg.east_top,:);
            bez.Ce4(reg.se_top,:) = newC4(reg.se_top,:);
            
            if norm(bez.Ce6(reg.ne_low,:))>0
                saveCell{6}=[saveCell{6}, 9];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.ne_low)];
            end
            
            if norm(bez.Ce6(reg.east_low,:))>0
                saveCell{6}=[saveCell{6}, 6];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.east_low)];
            end
            
            if norm(bez.Ce8(reg.east_low,:))>0
                saveCell{8}=[saveCell{8}, 6];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.east_low)];
            end
            
            if norm(bez.Ce8(reg.se_low,:))>0
                saveCell{8}=[saveCell{8}, 3];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.se_low)];
            end
            
            bez.Ce6(reg.ne_low,:) = newC6(reg.ne_low,:);
            bez.Ce6(reg.east_low,:) = newC6(reg.east_low,:);
            bez.Ce8(reg.east_low,:) = newC8(reg.east_low,:);
            bez.Ce8(reg.se_low,:) = newC8(reg.se_low,:);
            
            if norm(bez.Ce6(reg.ne_mid,:))>0
                saveCell{6}=[saveCell{6}, 18];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.ne_mid)];
            end
            
            if norm(bez.Ce6(reg.east_mid,:))>0
                saveCell{6}=[saveCell{6}, 15];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.east_mid)];
            end
            
            if norm(bez.Ce8(reg.east_mid,:))>0
                saveCell{8}=[saveCell{8}, 15];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.east_mid)];
            end
            
            if norm(bez.Ce8(reg.se_mid,:))>0
                saveCell{8}=[saveCell{8}, 12];
                savedIndex{6} = [savedIndex{8}, nodes.in8(reg.se_mid)];
            end
            
            bez.Ce6(reg.ne_mid,:) = newC6(reg.ne_mid,:);
            bez.Ce6(reg.east_mid,:) = newC6(reg.east_mid,:);
            bez.Ce8(reg.east_mid,:) = newC8(reg.east_mid,:);
            bez.Ce8(reg.se_mid,:) = newC8(reg.se_mid,:);
        case 3
            if norm(bez.Ce3(reg.ne_mid,:))>0
                saveCell{3}=[saveCell{3}, 18];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.ne_mid)];
            end
            if norm(bez.Ce3(reg.north_mid,:))>0
                saveCell{3}=[saveCell{3}, 17];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.north_mid)];
            end
            if norm(bez.Ce4(reg.north_mid,:))>0
                saveCell{4}=[saveCell{4}, 17];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.north_mid)];
            end
            if norm(bez.Ce4(reg.nw_mid,:))>0
                saveCell{4}=[saveCell{4}, 16];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.nw_mid)];
            end
            
            bez.Ce3(reg.ne_mid,:) = newC3(reg.ne_mid,:);
            bez.Ce3(reg.north_mid,:) = newC3(reg.north_mid,:);
            bez.Ce4(reg.north_mid,:) = newC4(reg.north_mid,:);
            bez.Ce4(reg.nw_mid,:) = newC4(reg.nw_mid,:);
            
            if norm(bez.Ce3(reg.ne_top,:))>0
                saveCell{3}=[saveCell{3}, 27];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.ne_top)];
            end
            if norm(bez.Ce3(reg.north_top,:))>0
                saveCell{3}=[saveCell{3}, 26];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.north_top)];
            end
            if norm(bez.Ce4(reg.north_top,:))>0
                saveCell{4}=[saveCell{4}, 26];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.north_top)];
            end
            if norm(bez.Ce4(reg.nw_top,:))>0
                saveCell{4}=[saveCell{4}, 25];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.nw_top)];
            end            
            
            bez.Ce3(reg.ne_top,:) = newC3(reg.ne_top,:);
            bez.Ce3(reg.north_top,:) = newC3(reg.north_top,:);
            bez.Ce4(reg.north_top,:) = newC4(reg.north_top,:);
            bez.Ce4(reg.nw_top,:) = newC4(reg.nw_top,:);
            
            if norm(bez.Ce7(reg.ne_low,:))>0
                saveCell{7}=[saveCell{7}, 9];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.ne_low)];
            end
            if norm(bez.Ce7(reg.north_low,:))>0
                saveCell{7}=[saveCell{7}, 8];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.north_low)];
            end
            if norm(bez.Ce8(reg.north_low,:))>0
                saveCell{8}=[saveCell{8}, 8];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.north_low)];
            end
            if norm(bez.Ce8(reg.nw_low,:))>0
                saveCell{8}=[saveCell{8}, 7];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.nw_low)];
            end
            
            bez.Ce7(reg.ne_low,:) = newC7(reg.ne_low,:);
            bez.Ce7(reg.north_low,:) = newC7(reg.north_low,:);
            bez.Ce8(reg.north_low,:) = newC8(reg.north_low,:);
            bez.Ce8(reg.nw_low,:) = newC8(reg.nw_low,:);
            
            if norm(bez.Ce7(reg.ne_mid,:))>0
                saveCell{7}=[saveCell{7}, 18];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.ne_mid)];
            end
            if norm(bez.Ce7(reg.north_mid,:))>0
                saveCell{7}=[saveCell{7}, 17];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.north_mid)];
            end
            if norm(bez.Ce8(reg.north_mid,:))>0
                saveCell{8}=[saveCell{8}, 17];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.north_mid)];
            end
            if norm(bez.Ce8(reg.nw_mid,:))>0
                saveCell{8}=[saveCell{8}, 16];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.nw_mid)];
            end
            
            bez.Ce7(reg.ne_mid,:) = newC7(reg.ne_mid,:);
            bez.Ce7(reg.north_mid,:) = newC7(reg.north_mid,:);
            bez.Ce8(reg.north_mid,:) = newC8(reg.north_mid,:);
            bez.Ce8(reg.nw_mid,:) = newC8(reg.nw_mid,:);
        case 4
            if norm(bez.Ce1(reg.nw_mid,:))>0
                saveCell{1}=[saveCell{1}, 16];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.nw_mid)];
            end
            
            if norm(bez.Ce1(reg.west_mid,:))>0
                saveCell{1}=[saveCell{1}, 13];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.west_mid)];
            end
            
            if norm(bez.Ce3(reg.west_mid,:))>0
                saveCell{3}=[saveCell{3}, 13];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.west_mid)];
            end
            
            if norm(bez.Ce3(reg.sw_mid,:))>0
                saveCell{3}=[saveCell{3}, 10];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.sw_mid)];
            end
            
            
            bez.Ce1(reg.nw_mid,:) = newC1(reg.nw_mid,:);
            bez.Ce1(reg.west_mid,:) = newC1(reg.west_mid,:);
            bez.Ce3(reg.west_mid,:) = newC3(reg.west_mid,:);
            bez.Ce3(reg.sw_mid,:) = newC3(reg.sw_mid,:);
            
            if norm(bez.Ce1(reg.nw_top,:))>0
                saveCell{1}=[saveCell{1}, 25];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.nw_top)];
            end
            
            if norm(bez.Ce1(reg.west_top,:))>0
                saveCell{1}=[saveCell{1}, 22];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.west_top)];
            end
            
            if norm(bez.Ce3(reg.west_top,:))>0
                saveCell{3}=[saveCell{3}, 22];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.west_top)];
            end
            
            if norm(bez.Ce3(reg.sw_top,:))>0
                saveCell{3}=[saveCell{3}, 19];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.sw_top)];
            end
            
            bez.Ce1(reg.nw_top,:) = newC1(reg.nw_top,:);
            bez.Ce1(reg.west_top,:) = newC1(reg.west_top,:);
            bez.Ce3(reg.west_top,:) = newC3(reg.west_top,:);
            bez.Ce3(reg.sw_top,:) = newC3(reg.sw_top,:);
            
            if norm(bez.Ce5(reg.nw_low,:))>0
                saveCell{5}=[saveCell{5}, 7];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.nw_low)];
            end
            
            if norm(bez.Ce5(reg.west_low,:))>0
                saveCell{5}=[saveCell{5}, 4];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.west_low)];
            end
            
            if norm(bez.Ce7(reg.west_low,:))>0
                saveCell{7}=[saveCell{7}, 4];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.west_low)];
            end
            
            if norm(bez.Ce7(reg.sw_low,:))>0
                saveCell{7}=[saveCell{7}, 1];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.sw_low)];
            end
            
            bez.Ce5(reg.nw_low,:) = newC5(reg.nw_low,:);
            bez.Ce5(reg.west_low,:) = newC5(reg.west_low,:);
            bez.Ce7(reg.west_low,:) = newC7(reg.west_low,:);
            bez.Ce7(reg.sw_low,:) = newC7(reg.sw_low,:);
            
            if norm(bez.Ce5(reg.nw_mid,:))>0
                saveCell{5}=[saveCell{5}, 16];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.nw_mid)];
            end
            
            if norm(bez.Ce5(reg.west_mid,:))>0
                saveCell{5}=[saveCell{5}, 13];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.west_mid)];
            end
            
            if norm(bez.Ce7(reg.west_mid,:))>0
                saveCell{7}=[saveCell{7}, 13];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.west_mid)];
            end
            
            if norm(bez.Ce7(reg.sw_mid,:))>0
                saveCell{7}=[saveCell{7}, 10];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.sw_mid)];
            end
            
            bez.Ce5(reg.nw_mid,:) = newC5(reg.nw_mid,:);
            bez.Ce5(reg.west_mid,:) = newC5(reg.west_mid,:);
            bez.Ce7(reg.west_mid,:) = newC7(reg.west_mid,:);
            bez.Ce7(reg.sw_mid,:) = newC7(reg.sw_mid,:);
        case 5
            if norm(bez.Ce1(reg.center_low,:))>0
                saveCell{1}=[saveCell{1}, 5];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.center_low)];
            end
            if norm(bez.Ce1(reg.east_low,:))>0
                saveCell{1}=[saveCell{1}, 6];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.east_low)];
            end
            if norm(bez.Ce1(reg.north_low,:))>0
                saveCell{1}=[saveCell{1}, 8];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.north_low)];
            end
            if norm(bez.Ce1(reg.ne_low,:))>0
                saveCell{1}=[saveCell{1}, 9];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.ne_low)];
            end
            bez.Ce1(reg.center_low,:) = newC1(reg.center_low,:);
            bez.Ce1(reg.east_low,:) = newC1(reg.east_low,:);
            bez.Ce1(reg.north_low,:) = newC1(reg.north_low,:);
            bez.Ce1(reg.ne_low,:) = newC1(reg.ne_low,:);
            
            if norm(bez.Ce2(reg.center_low,:))>0
                saveCell{2}=[saveCell{2}, 5];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.center_low)];
            end
            
            if norm(bez.Ce2(reg.west_low,:))>0
                saveCell{2}=[saveCell{2}, 4];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.west_low)];
            end
            
            if norm(bez.Ce2(reg.north_low,:))>0
                saveCell{2}=[saveCell{2}, 8];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.north_low)];
            end
            
            if norm(bez.Ce2(reg.nw_low,:))>0
                saveCell{2}=[saveCell{2}, 7];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.nw_low)];
            end
            
            bez.Ce2(reg.center_low,:) = newC2(reg.center_low,:);
            bez.Ce2(reg.west_low,:) = newC2(reg.west_low,:);
            bez.Ce2(reg.north_low,:) = newC2(reg.north_low,:);
            bez.Ce2(reg.nw_low,:) = newC2(reg.nw_low,:);
            
            if norm(bez.Ce3(reg.center_low,:))>0
                saveCell{3}=[saveCell{3}, 5];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.center_low)];
            end
            
            if norm(bez.Ce3(reg.south_low,:))>0
                saveCell{3}=[saveCell{3}, 2];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.south_low)];
            end
            
            if norm(bez.Ce3(reg.east_low,:))>0
                saveCell{3}=[saveCell{3}, 6];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.east_low)];
            end
            
            if norm(bez.Ce3(reg.se_low,:))>0
                saveCell{3}=[saveCell{3}, 3];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.se_low)];
            end
            
            bez.Ce3(reg.center_low,:) = newC3(reg.center_low,:);
            bez.Ce3(reg.south_low,:) = newC3(reg.south_low,:);
            bez.Ce3(reg.east_low,:) = newC3(reg.east_low,:);
            bez.Ce3(reg.se_low,:) = newC3(reg.se_low,:);
            
            if norm(bez.Ce4(reg.center_low,:))>0
                saveCell{4}=[saveCell{4}, 5];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.center_low)];
            end
            
            if norm(bez.Ce4(reg.south_low,:))>0
                saveCell{4}=[saveCell{4}, 2];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.south_low)];
            end
            
            if norm(bez.Ce4(reg.west_low,:))>0
                saveCell{4}=[saveCell{4}, 4];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.west_low)];
            end
            
            if norm(bez.Ce4(reg.sw_low,:))>0
                saveCell{4}=[saveCell{4}, 1];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.sw_low)];
            end
            
            bez.Ce4(reg.center_low,:) = newC4(reg.center_low,:);
            bez.Ce4(reg.south_low,:) = newC4(reg.south_low,:);
            bez.Ce4(reg.west_low,:) = newC4(reg.west_low,:);
            bez.Ce4(reg.sw_low,:) = newC4(reg.sw_low,:);
        case 6            
            if norm(bez.Ce5(reg.center_top,:))>0
                saveCell{5}=[saveCell{5}, 23];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.center_top)];
            end
            
            if norm(bez.Ce5(reg.east_top,:))>0
                saveCell{5}=[saveCell{5}, 24];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.east_top)];
            end
            
            if norm(bez.Ce5(reg.north_top,:))>0
                saveCell{5}=[saveCell{5}, 26];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.north_top)];
            end            
            if norm(bez.Ce5(reg.ne_top,:))>0
                saveCell{5}=[saveCell{5}, 27];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.ne_top)];
            end
            
            bez.Ce5(reg.center_top,:) = newC5(reg.center_top,:);
            bez.Ce5(reg.east_top,:) = newC5(reg.east_top,:);
            bez.Ce5(reg.north_top,:) = newC5(reg.north_top,:);
            bez.Ce5(reg.ne_top,:) = newC5(reg.ne_top,:);
            
            if norm(bez.Ce6(reg.center_top,:))>0
                saveCell{6}=[saveCell{6}, 23];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.center_top)];
            end
            
            if norm(bez.Ce6(reg.west_top,:))>0
                saveCell{6}=[saveCell{6}, 22];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.west_top)];
            end
            
            if norm(bez.Ce6(reg.north_top,:))>0
                saveCell{6}=[saveCell{6}, 26];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.north_top)];
            end
            
            if norm(bez.Ce6(reg.nw_top,:))>0
                saveCell{6}=[saveCell{6}, 25];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.nw_top)];
            end
            
            bez.Ce6(reg.center_top,:) = newC6(reg.center_top,:);
            bez.Ce6(reg.west_top,:) = newC6(reg.west_top,:);
            bez.Ce6(reg.north_top,:) = newC6(reg.north_top,:);
            bez.Ce6(reg.nw_top,:) = newC6(reg.nw_top,:);
            
            if norm(bez.Ce7(reg.center_top,:))>0
                saveCell{7}=[saveCell{7}, 23];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.center_top)];
            end
            
            if norm(bez.Ce7(reg.south_top,:))>0
                saveCell{7}=[saveCell{7}, 20];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.south_top)];
            end
            
            if norm(bez.Ce7(reg.east_top,:))>0
                saveCell{7}=[saveCell{7}, 24];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.east_top)];
            end
            
            if norm(bez.Ce7(reg.se_top,:))>0
                saveCell{7}=[saveCell{7}, 21];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.se_top)];
            end
            
            bez.Ce7(reg.center_top,:) = newC7(reg.center_top,:);
            bez.Ce7(reg.south_top,:) = newC7(reg.south_top,:);
            bez.Ce7(reg.east_top,:) = newC7(reg.east_top,:);
            bez.Ce7(reg.se_top,:) = newC7(reg.se_top,:);
            
            if norm(bez.Ce8(reg.center_top,:))>0
                saveCell{8}=[saveCell{8}, 23];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.center_top)];
            end
            
            if norm(bez.Ce8(reg.south_top,:))>0
                saveCell{8}=[saveCell{8}, 20];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.south_top)];
            end
            
            if norm(bez.Ce8(reg.west_top,:))>0
                saveCell{8}=[saveCell{8}, 22];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.west_top)];
            end
            
            if norm(bez.Ce8(reg.sw_top,:))>0
                saveCell{8}=[saveCell{8}, 19];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.sw_top)];
            end
            
            bez.Ce8(reg.center_top,:) = newC8(reg.center_top,:);
            bez.Ce8(reg.south_top,:) = newC8(reg.south_top,:);
            bez.Ce8(reg.west_top,:) = newC8(reg.west_top,:);
            bez.Ce8(reg.sw_top,:) = newC8(reg.sw_top,:);
        case 7
            if norm(bez.Ce1(reg.center_mid,:))>0
                saveCell{1}=[saveCell{1}, 14];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.center_mid)];
            end
            if norm(bez.Ce1(reg.east_mid,:))>0
                saveCell{1}=[saveCell{1}, 15];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.east_mid)];
            end
            if norm(bez.Ce1(reg.north_mid,:))>0
                saveCell{1}=[saveCell{1}, 17];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.north_mid)];
            end
            if norm(bez.Ce1(reg.ne_mid,:))>0
                saveCell{1}=[saveCell{1}, 18];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.ne_mid)];
            end
            if norm(bez.Ce2(reg.center_mid,:))>0
                saveCell{2}=[saveCell{2}, 14];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.center_mid)];
            end
            if norm(bez.Ce2(reg.west_mid,:))>0
                saveCell{2}=[saveCell{2}, 13];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.west_mid)];
            end
            if norm(bez.Ce2(reg.north_mid,:))>0
                saveCell{2}=[saveCell{2}, 17];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.north_mid)];
            end
            if norm(bez.Ce2(reg.nw_mid,:))>0
                saveCell{2}=[saveCell{2}, 16];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.nw_mid)];
            end
            if norm(bez.Ce3(reg.center_mid,:))>0
                saveCell{3}=[saveCell{3}, 14];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.center_mid)];
            end
            if norm(bez.Ce3(reg.south_mid,:))>0
                saveCell{3}=[saveCell{3}, 11];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.south_mid)];
            end
            if norm(bez.Ce3(reg.east_mid,:))>0
                saveCell{3}=[saveCell{3}, 15];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.east_mid)];
            end
            if norm(bez.Ce3(reg.se_mid,:))>0
                saveCell{3}=[saveCell{3}, 12];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.se_mid)];
            end
            if norm(bez.Ce4(reg.center_mid,:))>0
                saveCell{4}=[saveCell{4}, 14];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.center_mid)];
            end
            if norm(bez.Ce4(reg.south_mid,:))>0
                saveCell{4}=[saveCell{4}, 11];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.south_mid)];
            end
            if norm(bez.Ce4(reg.west_mid,:))>0
                saveCell{4}=[saveCell{4}, 13];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.west_mid)];
            end
            if norm(bez.Ce4(reg.sw_mid,:))>0
                saveCell{4}=[saveCell{4}, 10];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.sw_mid)];
            end
            
            %%%%%
            if norm(bez.Ce1(reg.center_top,:))>0
                saveCell{1}=[saveCell{1}, 23];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.center_top)];
            end
            if norm(bez.Ce1(reg.east_top,:))>0
                saveCell{1}=[saveCell{1}, 24];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.east_top)];
            end
            if norm(bez.Ce1(reg.north_top,:))>0
                saveCell{1}=[saveCell{1}, 26];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.north_top)];
            end
            if norm(bez.Ce1(reg.ne_top,:))>0
                saveCell{1}=[saveCell{1}, 27];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.ne_top)];
            end
            if norm(bez.Ce2(reg.center_top,:))>0
                saveCell{2}=[saveCell{2}, 23];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.center_top)];
            end
            if norm(bez.Ce2(reg.west_top,:))>0
                saveCell{2}=[saveCell{2}, 22];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.west_top)];
            end
            if norm(bez.Ce2(reg.north_top,:))>0
                saveCell{2}=[saveCell{2}, 26];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.north_top)];
            end
            if norm(bez.Ce2(reg.nw_top,:))>0
                saveCell{2}=[saveCell{2}, 25];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.nw_top)];
            end
            if norm(bez.Ce3(reg.center_top,:))>0
                saveCell{3}=[saveCell{3}, 23];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.center_top)];
            end
            if norm(bez.Ce3(reg.south_top,:))>0
                saveCell{3}=[saveCell{3}, 20];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.south_top)];
            end
            if norm(bez.Ce3(reg.east_top,:))>0
                saveCell{3}=[saveCell{3}, 24];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.east_top)];
            end
            if norm(bez.Ce3(reg.se_top,:))>0
                saveCell{3}=[saveCell{3}, 21];
                savedIndex{3} = [savedIndex{3}, nodes.in3(reg.se_top)];
            end
            if norm(bez.Ce4(reg.center_top,:))>0
                saveCell{4}=[saveCell{4}, 23];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.center_top)];
            end
            if norm(bez.Ce4(reg.south_top,:))>0
                saveCell{4}=[saveCell{4}, 20];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.south_top)];
            end
            if norm(bez.Ce4(reg.west_top,:))>0
                saveCell{4}=[saveCell{4}, 22];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.west_top)];
            end
            if norm(bez.Ce4(reg.sw_top,:))>0
                saveCell{4}=[saveCell{4}, 19];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.sw_top)];
            end
            
            %%
            if norm(bez.Ce5(reg.center_low,:))>0
                saveCell{5}=[saveCell{5}, 5];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.center_low)];
            end
            if norm(bez.Ce5(reg.east_low,:))>0
                saveCell{5}=[saveCell{5}, 6];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.east_low)];
            end
            if norm(bez.Ce5(reg.north_low,:))>0
                saveCell{5}=[saveCell{5}, 8];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.north_low)];
            end
            if norm(bez.Ce5(reg.ne_low,:))>0
                saveCell{5}=[saveCell{5}, 9];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.ne_low)];
            end
            if norm(bez.Ce6(reg.center_low,:))>0
                saveCell{6}=[saveCell{6}, 5];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.center_low)];
            end
            if norm(bez.Ce6(reg.west_low,:))>0
                saveCell{6}=[saveCell{6}, 4];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.west_low)];
            end
            if norm(bez.Ce6(reg.north_low,:))>0
                saveCell{6}=[saveCell{6}, 8];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.north_low)];
            end
            if norm(bez.Ce6(reg.nw_low,:))>0
                saveCell{6}=[saveCell{6}, 7];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.nw_low)];
            end
            if norm(bez.Ce7(reg.center_low,:))>0
                saveCell{7}=[saveCell{7}, 5];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.center_low)];
            end
            if norm(bez.Ce7(reg.south_low,:))>0
                saveCell{7}=[saveCell{7}, 2];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.south_low)];
            end
            if norm(bez.Ce7(reg.east_low,:))>0
                saveCell{7}=[saveCell{7}, 6];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.east_low)];
            end
            if norm(bez.Ce7(reg.se_low,:))>0
                saveCell{7}=[saveCell{7}, 3];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.se_low)];
            end
            if norm(bez.Ce8(reg.center_low,:))>0
                saveCell{8}=[saveCell{8}, 5];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.center_low)];
            end
            if norm(bez.Ce8(reg.south_low,:))>0
                saveCell{8}=[saveCell{8}, 2];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.south_low)];
            end
            if norm(bez.Ce8(reg.west_low,:))>0
                saveCell{8}=[saveCell{8}, 4];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.west_low)];
            end
            if norm(bez.Ce8(reg.sw_low,:))>0
                saveCell{8}=[saveCell{8}, 1];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.sw_low)];
            end
            
            %%
            
            if norm(bez.Ce5(reg.center_mid,:))>0
                saveCell{5}=[saveCell{5}, 14];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.center_mid)];
            end
            if norm(bez.Ce5(reg.east_mid,:))>0
                saveCell{5}=[saveCell{5}, 15];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.east_mid)];
            end
            if norm(bez.Ce5(reg.north_mid,:))>0
                saveCell{5}=[saveCell{5}, 17];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.north_mid)];
            end
            if norm(bez.Ce5(reg.ne_mid,:))>0
                saveCell{5}=[saveCell{5}, 18];
                savedIndex{5} = [savedIndex{5}, nodes.in5(reg.ne_mid)];
            end
            if norm(bez.Ce6(reg.center_mid,:))>0
                saveCell{6}=[saveCell{6}, 14];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.center_mid)];
            end
            if norm(bez.Ce6(reg.west_mid,:))>0
                saveCell{6}=[saveCell{6}, 13];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.west_mid)];
            end
            if norm(bez.Ce6(reg.north_mid,:))>0
                saveCell{6}=[saveCell{6}, 17];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.north_mid)];
            end
            if norm(bez.Ce6(reg.nw_mid,:))>0
                saveCell{6}=[saveCell{6}, 16];
                savedIndex{6} = [savedIndex{6}, nodes.in6(reg.nw_mid)];
            end
            if norm(bez.Ce7(reg.center_mid,:))>0
                saveCell{7}=[saveCell{7}, 14];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.center_mid)];
            end
            if norm(bez.Ce7(reg.south_mid,:))>0
                saveCell{7}=[saveCell{7}, 11];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.south_mid)];
            end
            if norm(bez.Ce7(reg.east_mid,:))>0
                saveCell{7}=[saveCell{7}, 15];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.east_mid)];
            end
            if norm(bez.Ce7(reg.se_mid,:))>0
                saveCell{7}=[saveCell{7}, 12];
                savedIndex{7} = [savedIndex{7}, nodes.in7(reg.se_mid)];
            end
            if norm(bez.Ce8(reg.center_mid,:))>0
                saveCell{8}=[saveCell{8}, 14];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.center_mid)];
            end
            if norm(bez.Ce8(reg.south_mid,:))>0
                saveCell{8}=[saveCell{8}, 11];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.south_mid)];
            end
            if norm(bez.Ce8(reg.west_mid,:))>0
                saveCell{8}=[saveCell{8}, 13];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.west_mid)];
            end
            if norm(bez.Ce8(reg.sw_mid,:))>0
                saveCell{8}=[saveCell{8}, 10];
                savedIndex{8} = [savedIndex{8}, nodes.in8(reg.sw_mid)];
            end
            
            %%
            
            %mid_slice in elements 1,2,3,4
            bez.Ce1(reg.center_mid,:) = newC1(reg.center_mid,:);
            bez.Ce1(reg.east_mid,:) = newC1(reg.east_mid,:);
            bez.Ce1(reg.north_mid,:) = newC1(reg.north_mid,:);
            bez.Ce1(reg.ne_mid,:) = newC1(reg.ne_mid,:);
            bez.Ce2(reg.center_mid,:) = newC2(reg.center_mid,:);
            bez.Ce2(reg.west_mid,:) = newC2(reg.west_mid,:);
            bez.Ce2(reg.north_mid,:) = newC2(reg.north_mid,:);
            bez.Ce2(reg.nw_mid,:) = newC2(reg.nw_mid,:);
            bez.Ce3(reg.center_mid,:) = newC3(reg.center_mid,:);
            bez.Ce3(reg.south_mid,:) = newC3(reg.south_mid,:);
            bez.Ce3(reg.east_mid,:) = newC3(reg.east_mid,:);
            bez.Ce3(reg.se_mid,:) = newC3(reg.se_mid,:);
            bez.Ce4(reg.center_mid,:) = newC4(reg.center_mid,:);
            bez.Ce4(reg.south_mid,:) = newC4(reg.south_mid,:);
            bez.Ce4(reg.west_mid,:) = newC4(reg.west_mid,:);
            bez.Ce4(reg.sw_mid,:) = newC4(reg.sw_mid,:);
            %top_slice in elements 1,2,3,4
            bez.Ce1(reg.center_top,:) = newC1(reg.center_top,:);
            bez.Ce1(reg.east_top,:) = newC1(reg.east_top,:);
            bez.Ce1(reg.north_top,:) = newC1(reg.north_top,:);
            bez.Ce1(reg.ne_top,:) = newC1(reg.ne_top,:);
            bez.Ce2(reg.center_top,:) = newC2(reg.center_top,:);
            bez.Ce2(reg.west_top,:) = newC2(reg.west_top,:);
            bez.Ce2(reg.north_top,:) = newC2(reg.north_top,:);
            bez.Ce2(reg.nw_top,:) = newC2(reg.nw_top,:);
            bez.Ce3(reg.center_top,:) = newC3(reg.center_top,:);
            bez.Ce3(reg.south_top,:) = newC3(reg.south_top,:);
            bez.Ce3(reg.east_top,:) = newC3(reg.east_top,:);
            bez.Ce3(reg.se_top,:) = newC3(reg.se_top,:);
            bez.Ce4(reg.center_top,:) = newC4(reg.center_top,:);
            bez.Ce4(reg.south_top,:) = newC4(reg.south_top,:);
            bez.Ce4(reg.west_top,:) = newC4(reg.west_top,:);
            bez.Ce4(reg.sw_top,:) = newC4(reg.sw_top,:);
            %low_slice in elements 5,6,7,8
            bez.Ce5(reg.center_low,:) = newC5(reg.center_low,:);
            bez.Ce5(reg.east_low,:) = newC5(reg.east_low,:);
            bez.Ce5(reg.north_low,:) = newC5(reg.north_low,:);
            bez.Ce5(reg.ne_low,:) = newC5(reg.ne_low,:);
            bez.Ce6(reg.center_low,:) = newC6(reg.center_low,:);
            bez.Ce6(reg.west_low,:) = newC6(reg.west_low,:);
            bez.Ce6(reg.north_low,:) = newC6(reg.north_low,:);
            bez.Ce6(reg.nw_low,:) = newC6(reg.nw_low,:);
            bez.Ce7(reg.center_low,:) = newC7(reg.center_low,:);
            bez.Ce7(reg.south_low,:) = newC7(reg.south_low,:);
            bez.Ce7(reg.east_low,:) = newC7(reg.east_low,:);
            bez.Ce7(reg.se_low,:) = newC7(reg.se_low,:);
            bez.Ce8(reg.center_low,:) = newC8(reg.center_low,:);
            bez.Ce8(reg.south_low,:) = newC8(reg.south_low,:);
            bez.Ce8(reg.west_low,:) = newC8(reg.west_low,:);
            bez.Ce8(reg.sw_low,:) = newC8(reg.sw_low,:);
            %mid_slice in elements 5,6,7,8
            bez.Ce5(reg.center_mid,:) = newC5(reg.center_mid,:);
            bez.Ce5(reg.east_mid,:) = newC5(reg.east_mid,:);
            bez.Ce5(reg.north_mid,:) = newC5(reg.north_mid,:);
            bez.Ce5(reg.ne_mid,:) = newC5(reg.ne_mid,:);
            bez.Ce6(reg.center_mid,:) = newC6(reg.center_mid,:);
            bez.Ce6(reg.west_mid,:) = newC6(reg.west_mid,:);
            bez.Ce6(reg.north_mid,:) = newC6(reg.north_mid,:);
            bez.Ce6(reg.nw_mid,:) = newC6(reg.nw_mid,:);
            bez.Ce7(reg.center_mid,:) = newC7(reg.center_mid,:);
            bez.Ce7(reg.south_mid,:) = newC7(reg.south_mid,:);
            bez.Ce7(reg.east_mid,:) = newC7(reg.east_mid,:);
            bez.Ce7(reg.se_mid,:) = newC7(reg.se_mid,:);
            bez.Ce8(reg.center_mid,:) = newC8(reg.center_mid,:);
            bez.Ce8(reg.south_mid,:) = newC8(reg.south_mid,:);
            bez.Ce8(reg.west_mid,:) = newC8(reg.west_mid,:);
            bez.Ce8(reg.sw_mid,:) = newC8(reg.sw_mid,:);
        case 8
            if norm(bez.Ce5(reg.west_top,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce5(reg.nw_top,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce7(reg.sw_top,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce7(reg.west_top,:))>0
                disp(newVert(j))
            end
            bez.Ce5(reg.west_top,:) = newC5(reg.west_top,:);
            bez.Ce5(reg.nw_top,:) = newC5(reg.nw_top,:);
            bez.Ce7(reg.sw_top,:) = newC7(reg.sw_top,:);
            bez.Ce7(reg.west_top,:) = newC7(reg.west_top,:);
        case 9
            if norm(bez.Ce1(reg.west_low,:))>0                
                saveCell{1}=[saveCell{1}, 4];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.west_low)];
            end
            if norm(bez.Ce1(reg.nw_low,:))>0                
                saveCell{1}=[saveCell{1}, 7];
                savedIndex{1} = [savedIndex{1}, nodes.in1(reg.nw_low)];
            end
            if norm(bez.Ce3(reg.sw_low,:))>0
               saveCell{3}=[saveCell{3}, 1];
               savedIndex{3} = [savedIndex{3}, nodes.in3(reg.sw_low)];
            end
            if norm(bez.Ce3(reg.west_low,:))>0
               saveCell{3}=[saveCell{3}, 4];
               savedIndex{3} = [savedIndex{3}, nodes.in3(reg.west_low)];
            end
            bez.Ce1(reg.west_low,:) = newC1(reg.west_low,:);
            bez.Ce1(reg.nw_low,:) = newC1(reg.nw_low,:);
            bez.Ce3(reg.sw_low,:) = newC3(reg.sw_low,:);
            bez.Ce3(reg.west_low,:) = newC3(reg.west_low,:);
        case 10
            if norm(bez.Ce6(reg.east_top,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce6(reg.ne_top,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce8(reg.se_top,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce8(reg.east_top,:))>0
                disp(newVert(j))
            end
            bez.Ce6(reg.east_top,:) = newC6(reg.east_top,:);
            bez.Ce6(reg.ne_top,:) = newC6(reg.ne_top,:);
            bez.Ce8(reg.se_top,:) = newC8(reg.se_top,:);
            bez.Ce8(reg.east_top,:) = newC8(reg.east_top,:);
        case 11
            if norm(bez.Ce2(reg.east_low,:))>0
                saveCell{2}=[saveCell{2}, 6];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.east_low)];
            end
            if norm(bez.Ce2(reg.ne_low,:))>0                
                saveCell{2}=[saveCell{2}, 9];
                savedIndex{2} = [savedIndex{2}, nodes.in2(reg.ne_low)];
            end
            if norm(bez.Ce4(reg.se_low,:))>0
                saveCell{4}=[saveCell{4}, 3];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.se_low)];
            end
            if norm(bez.Ce4(reg.east_low,:))>0
                saveCell{4}=[saveCell{4}, 6];
                savedIndex{4} = [savedIndex{4}, nodes.in4(reg.east_low)];
            end
            bez.Ce2(reg.east_low,:) = newC2(reg.east_low,:);
            bez.Ce2(reg.ne_low,:) = newC2(reg.ne_low,:);
            bez.Ce4(reg.se_low,:) = newC4(reg.se_low,:);
            bez.Ce4(reg.east_low,:) = newC4(reg.east_low,:);
        case 12
            if norm(bez.Ce5(reg.south_top,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce5(reg.se_top,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce6(reg.sw_top,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce6(reg.south_top,:))>0
                disp(newVert(j))
            end
            bez.Ce5(reg.south_top,:) = newC5(reg.south_top,:);
            bez.Ce5(reg.se_top,:) = newC5(reg.se_top,:);
            bez.Ce6(reg.sw_top,:) = newC6(reg.sw_top,:);
            bez.Ce6(reg.south_top,:) = newC6(reg.south_top,:);
        case 13
            if norm(bez.Ce1(reg.south_low,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce1(reg.se_low,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce2(reg.sw_low,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce2(reg.south_low,:))>0
                disp(newVert(j))
            end
            bez.Ce1(reg.south_low,:) = newC1(reg.south_low,:);
            bez.Ce1(reg.se_low,:) = newC1(reg.se_low,:);
            bez.Ce2(reg.sw_low,:) = newC2(reg.sw_low,:);
            bez.Ce2(reg.south_low,:) = newC2(reg.south_low,:);
        case 14
            if norm(bez.Ce7(reg.north_top,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce7(reg.ne_top,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce8(reg.nw_top,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce8(reg.north_top,:))>0
                disp(newVert(j))
            end
            bez.Ce7(reg.north_top,:) = newC7(reg.north_top,:);
            bez.Ce7(reg.ne_top,:) = newC7(reg.ne_top,:);
            bez.Ce8(reg.nw_top,:) = newC8(reg.nw_top,:);
            bez.Ce8(reg.north_top,:) = newC8(reg.north_top,:);
        case 15
            if norm(bez.Ce3(reg.north_low,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce3(reg.ne_low,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce4(reg.nw_low,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce4(reg.north_low,:))>0
                disp(newVert(j))
            end
            bez.Ce3(reg.north_low,:) = newC3(reg.north_low,:);
            bez.Ce3(reg.ne_low,:) = newC3(reg.ne_low,:);
            bez.Ce4(reg.nw_low,:) = newC4(reg.nw_low,:);
            bez.Ce4(reg.north_low,:) = newC4(reg.north_low,:);
        case 16
            if norm(bez.Ce1(reg.sw_mid,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce1(reg.sw_top,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce5(reg.sw_low,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce5(reg.sw_mid,:))>0
                disp(newVert(j))
            end
            bez.Ce1(reg.sw_mid,:) = newC1(reg.sw_mid,:);
            bez.Ce1(reg.sw_top,:) = newC1(reg.sw_top,:);
            bez.Ce5(reg.sw_low,:) = newC5(reg.sw_low,:);
            bez.Ce5(reg.sw_mid,:) = newC5(reg.sw_mid,:);
        case 17
            if norm(bez.Ce2(reg.se_mid,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce2(reg.se_top,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce6(reg.se_low,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce6(reg.se_mid,:))>0
                disp(newVert(j))
            end
            bez.Ce2(reg.se_mid,:) = newC2(reg.se_mid,:);
            bez.Ce2(reg.se_top,:) = newC2(reg.se_top,:);
            bez.Ce6(reg.se_low,:) = newC6(reg.se_low,:);
            bez.Ce6(reg.se_mid,:) = newC6(reg.se_mid,:);
        case 18
            if norm(bez.Ce3(reg.nw_mid,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce3(reg.nw_top,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce7(reg.nw_low,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce7(reg.nw_mid,:))>0
                disp(newVert(j))
            end
            bez.Ce3(reg.nw_mid,:) = newC3(reg.nw_mid,:);
            bez.Ce3(reg.nw_top,:) = newC3(reg.nw_top,:);
            bez.Ce7(reg.nw_low,:) = newC7(reg.nw_low,:);
            bez.Ce7(reg.nw_mid,:) = newC7(reg.nw_mid,:);
        case 19
            if norm(bez.Ce4(reg.ne_mid,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce4(reg.ne_top,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce8(reg.ne_low,:))>0
                disp(newVert(j))
            end
            if norm(bez.Ce8(reg.ne_mid,:))>0
                disp(newVert(j))
            end
            bez.Ce4(reg.ne_mid,:) = newC4(reg.ne_mid,:);
            bez.Ce4(reg.ne_top,:) = newC4(reg.ne_top,:);
            bez.Ce8(reg.ne_low,:) = newC8(reg.ne_low,:);
            bez.Ce8(reg.ne_mid,:) = newC8(reg.ne_mid,:);
    end
end

%over-write the element_nod indices
nodes.in1 = newBasisSet{1};
nodes.in2 = newBasisSet{2};
nodes.in3 = newBasisSet{3};
nodes.in4 = newBasisSet{4};
nodes.in5 = newBasisSet{5};
nodes.in6 = newBasisSet{6};
nodes.in7 = newBasisSet{7};
nodes.in8 = newBasisSet{8};

%add the saved basis functions
Ce_extra = cell(1,8);
for save_index = 1:8
    for region_index = 1:length(saveCell{save_index})
        switch saveCell{save_index}(region_index)
            case 1  %sw_low corner
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.sw_low,:)];
            case 2   %south_low
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.south_low,:)];
            case 3  %se_low corner
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.se_low,:)];
            case 4  %west_low
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.west_low,:)];
            case 5  %center_low
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.center_low,:)];
            case 6  %east_low
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.east_low,:)];
            case 7  %nw_low corner
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.nw_low,:)];
            case 8  %north_low
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.north_low,:)];
            case 9  %ne_low corner
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.ne_low,:)];
            case 10  %sw_mid corner
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.sw_mid,:)];
            case 11   %south_mid
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.south_mid,:)];
            case 12  %se_mid corner
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.se_mid,:)];
            case 13  %west_mid
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.west_mid,:)];
            case 14  %center_mid
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.center_mid,:)];
            case 15  %east_mid
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.east_mid,:)];
            case 16  %nw_mid corner
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.nw_mid,:)];
            case 17  %north_mid
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.north_mid,:)];
            case 18  %ne_mid corner
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.ne_mid,:)];
            case 19  %sw_top corner
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.sw_top,:)];
            case 20   %south_top
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.south_top,:)];
            case 21  %se_top corner
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.se_top,:)];
            case 22  %west_top
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.west_top,:)];
            case 23  %center_top
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.center_top,:)];
            case 24  %east_top
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.east_top,:)];
            case 25  %nw_top corner
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.nw_top,:)];
            case 26  %north_top
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.north_top,:)];
            case 27  %ne_top corner
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(reg.ne_top,:)];
        end
    end
end


%return the outputs
bez.Ce1 = [bez.Ce1; Ce_extra{1}];
bez.Ce2 = [bez.Ce2; Ce_extra{2}];
bez.Ce3 = [bez.Ce3; Ce_extra{3}];
bez.Ce4 = [bez.Ce4; Ce_extra{4}];
bez.Ce5 = [bez.Ce5; Ce_extra{5}];
bez.Ce6 = [bez.Ce6; Ce_extra{6}];
bez.Ce7 = [bez.Ce7; Ce_extra{7}];
bez.Ce8 = [bez.Ce8; Ce_extra{8}];


nodes.in1 = [nodes.in1, savedIndex{1}];
nodes.in2 = [nodes.in2, savedIndex{2}];
nodes.in3 = [nodes.in3, savedIndex{3}];
nodes.in4 = [nodes.in4, savedIndex{4}];
nodes.in5 = [nodes.in5, savedIndex{5}];
nodes.in6 = [nodes.in6, savedIndex{6}];
nodes.in7 = [nodes.in7, savedIndex{7}];
nodes.in8 = [nodes.in8, savedIndex{8}];

% size(bez.Ce2)
% size(nodes.in2)
% 
% saveCell{2}
% pause