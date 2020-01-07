function [ row1, row2 ] = deCasteljau1d( row )
% Split row array into two arrays of the same size using deCasteljau
% algorithm ( http://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm )

num_ent = size(row,2);
row2 = zeros(1, num_ent);

temp = zeros(num_ent, num_ent);

temp(1,:) = row;
for j=1:num_ent-1
    temp(j+1,1:end-1) = (temp(j,1:end-1)+temp(j,2:end))/2;
end
row1 = temp(:,1)';
for j=1:num_ent
    row2(j) = temp(end+1-j,j);
end
