function [ temp ] = deCasteljau3DTP( cur_row, p, q, r )
temp = zeros(2*(q+1), 2*(p+1), 2*(r+1));
cur_row_cube = permute(reshape(cur_row,p+1,q+1,r+1),[2,1,3]);
%do the 1d deCasteljau algorithm in the row direction for each slice
for j=1:q+1
    for k=1:r+1
        [temp1, temp2] = deCasteljau1d(cur_row_cube(j,:,k));
        temp(j,:,k) = [temp1, temp2];
    end
end
%do the 1d deCasteljau algorithm in the column direction for each slice
for j=1:2*(p+1)
    for k=1:r+1
        [temp1, temp2] = deCasteljau1d(temp(1:(q+1),j,k)');
        temp(:,j,k) = [temp1, temp2];
    end
end

%do the 1d deCasteljau algorithm in the slice direction for each row
%and column
for j=1:2*(q+1)
    for k=1:2*(p+1)
        [temp1, temp2] = deCasteljau1d(squeeze(temp(j,k,1:(r+1)))');
        temp(j,k,:) = [temp1, temp2];
    end
end


