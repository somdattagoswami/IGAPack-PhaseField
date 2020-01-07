function [ Ce1, Ce2, Ce3, Ce4, in1, in2, in3, in4, dimBasis, newBasisIndex ] = deCasteljau2dai_coordinates( Ce, knotU1, knotU2, knotV1, knotV2, p, q, newVert, knotUl, knotUr, knotVd, knotVu, elmIn, dimBasis)
% Split element with Bezier extraction operator Ce into 4 elements with
% Bezier extraction operator Ce1, Ce2, Ce3, Ce4
% saveCell: cell array of dimension 4, indicating which entries are saved
%newVert: new Vertices
%automatic detection of extra basis vertices (saved entries)
%calculates element_nod indices for new elements
%encoding scheme:  
%  o-----3------o
%  |     |      |
%  4-----5------2
%  |     |      |
%  o-----1------o

alpha = floor((p-1)/2);
beta = floor((q-1)/2);

% alpha = 1;
% beta = 1;

num_rows = size(Ce,1);
num_cols = size(Ce,2);

Ce1 = zeros(num_rows, num_cols);
Ce2 = zeros(num_rows, num_cols);
Ce3 = zeros(num_rows, num_cols);
Ce4 = zeros(num_rows, num_cols);

knotUmid = (knotU1+knotU2)/2;
knotVmid = (knotV1+knotV2)/2;

newKnotU = [knotUl*ones(1,p+1), knotU1*ones(1,p-alpha), knotUmid*ones(1,p-alpha), knotU2*ones(1,p-alpha), knotUr*ones(1, p+1)];
newKnotV = [knotVd*ones(1,q+1), knotV1*ones(1,q-alpha), knotVmid*ones(1,q-alpha), knotV2*ones(1,q-alpha), knotVu*ones(1, q+1)];
[newC_u, ~] = bezierExtraction(newKnotU, p);
[newC_v, ~] = bezierExtraction(newKnotV, q);

%calculate the tensor product Bezier extraction operators on the four
%new elements
newC1 = kron(newC_v(:,:,2),newC_u(:,:,2));
newC2 = kron(newC_v(:,:,2),newC_u(:,:,3));
newC3 = kron(newC_v(:,:,3),newC_u(:,:,2));
newC4 = kron(newC_v(:,:,3),newC_u(:,:,3));     

Ce_save = cell(1,4);

%copy the basis function indices to the children elements
in1 = elmIn;
in2 = elmIn;
in3 = elmIn;
in4 = elmIn;

%do tensor product deCasteljau
for i=1:size(Ce,1)
    temp = zeros(2*(q+1), 2*(p+1));
    cur_row = Ce(i,:);
    
    cur_row_sq = reshape(cur_row,p+1,q+1)';       
    %do the 1d deCasteljau algorithm in the row direction
    for j=1:q+1
        [temp1, temp2] = deCasteljau1d(cur_row_sq(j,:));
        temp(j,:) = [temp1, temp2];
    end
    %do the 1d deCasteljau algorithm in the column direction
    for j=1:2*(p+1)
        [temp1, temp2] = deCasteljau1d(temp(1:(q+1),j)');
        temp(:,j) = [temp1'; temp2'];
    end          
    
    %zero out entries coresponding to new vertices
    for j=1:length(newVert)
        switch newVert(j)
            case 1
                temp(1:q-beta,(alpha+2):2*p-alpha+1) =0;
            case 2
                temp((beta+2):2*q-beta+1,p+alpha+2:end) = 0;
            case 3
                temp(q+beta+2:end,(alpha+2):2*p-alpha+1) = 0;
            case 4
                temp((beta+2):2*q-beta+1,1:p-alpha) = 0;
            case 5
                temp((beta+2):2*q-beta+1,(alpha+2):2*p-alpha+1) = 0;
        end
    end
    
%      zero out entries coresponding to new vertices
%     for j=1:length(newVert)
%         switch newVert(j)
%             case 1
%                 temp(1:q-1,3:2*p) =0;
%             case 2
%                 temp(3:2*q,end-p+1:end) = 0;
%             case 3
%                 temp(end-q+1:end,3:2*p) = 0;
%             case 4
%                 temp(3:2*q,1:p-1) = 0;
%             case 5
%                 temp(3:2*q,3:2*p) = 0;
%         end
%     end
         
    Ce1(i,:) = reshape(temp(1:q+1,1:p+1)',1,(p+1)*(q+1));
    Ce2(i,:) = reshape(temp(1:q+1,p+2:end)',1,(p+1)*(q+1));
    Ce3(i,:) = reshape(temp(q+2:end,1:p+1)',1,(p+1)*(q+1));
    Ce4(i,:) = reshape(temp(q+2:end,p+2:end)',1,(p+1)*(q+1));             
    
    Ce_save{1}(i,:) = Ce1(i,:);
    Ce_save{2}(i,:) = Ce2(i,:);
    Ce_save{3}(i,:) = Ce3(i,:);
    Ce_save{4}(i,:) = Ce4(i,:);    
                     
                   
end

%add in the new basis functions
%define corners indices
saveCell = cell(1, 4);
savedIndex = cell(1,4);


[ sw, south, se, west, center, east, nw, north, ne ] = getCornerIndices( p,q );
[newBasisSet, dimBasis, newBasisIndex ] = newBasisIndices_coordinates( newVert, elmIn, p, q, dimBasis );

for j=1:length(newVert)
    switch newVert(j)
        case 1
            if norm(Ce1(south,:))>0
                saveCell{1}=[saveCell{1},2]; 
                savedIndex{1} = [savedIndex{1}, in1(south)];
            end
            if norm(Ce1(se,:))>0
                saveCell{1}=[saveCell{1},3];
                savedIndex{1} = [savedIndex{1}, in1(se)];
            end
            if norm(Ce2(south,:))>0
                saveCell{2}=[saveCell{1},2]; 
                savedIndex{1} = [savedIndex{1}, in2(south)];
            end
            if norm(Ce2(sw,:))>0
                saveCell{2}=[saveCell{2},1];
                savedIndex{2} = [savedIndex{2}, in2(sw)];
            end
            Ce1(se,:) = newC1(se,:);
            Ce1(south,:) = newC1(south,:);
            Ce2(south,:) = newC2(south,:);
            Ce2(sw,:) = newC2(sw,:);        
            
        case 2
            if norm(Ce2(ne,:))>0
                saveCell{2}=[saveCell{2},9];            
                savedIndex{2} = [savedIndex{2}, in2(ne)];
            end
            if norm(Ce2(east,:))>0
                saveCell{2} = [saveCell{2},6];
                savedIndex{2} = [savedIndex{2}, in2(east)];
            end
            if norm(Ce4(east,:))>0
                saveCell{4} = [saveCell{4},6];
                savedIndex{4} = [savedIndex{4}, in4(east)];
            end
            if norm(Ce4(se,:))>0
                saveCell{4}=[saveCell{4},3];
                savedIndex{4} = [savedIndex{4}, in4(se)];
            end
            Ce2(ne,:) = newC2(ne,:);
            Ce2(east,:) = newC2(east,:);
            Ce4(se,:) = newC4(se,:);
            Ce4(east,:) = newC4(east,:);
            
            
        case 3
            if norm(Ce3(ne,:))>0
                saveCell{3}=[saveCell{3},9];            
                savedIndex{3} = [savedIndex{3}, in3(ne)];
            end
            if norm(Ce3(north,:))>0
                saveCell{3}=[saveCell{3},8];
                savedIndex{3}=[savedIndex{3}, in3(north)];
            end
            if norm(Ce4(north,:))>0
                saveCell{4}=[saveCell{4},8];
                savedIndex{4}=[savedIndex{4}, in4(north)];
            end
            if norm(Ce4(nw,:))>0
                saveCell{4}=[saveCell{4},7];
                savedIndex{4} = [savedIndex{4}, in4(nw)];
            end
            Ce3(ne,:) = newC3(ne,:);
            Ce3(north,:) = newC3(north,:);
            Ce4(nw,:) = newC4(nw,:);      
            Ce4(north,:) = newC4(north,:);
            
            
        case 4
            if norm(Ce1(nw,:))>0
                saveCell{1}=[saveCell{1},7];    
                savedIndex{1} = [savedIndex{1}, in1(nw)];
            end
            if norm(Ce1(west,:))>0
                saveCell{1}=[saveCell{1},4];
                savedIndex{1}=[savedIndex{1}, in1(west)];
            end
            if norm(Ce3(sw,:))>0
                saveCell{3}=[saveCell{3},1];
                savedIndex{3} = [savedIndex{3}, in3(sw)];
            end
            if norm(Ce3(west,:))>0
                saveCell{3}=[saveCell{3},4];
                savedIndex{3} = [savedIndex{3}, in3(west)];
            end
            Ce1(nw,:) = newC1(nw,:);
            Ce1(west,:) = newC1(west,:);
            Ce3(sw,:) = newC3(sw,:);    
            Ce3(west,:) = newC3(west,:);
            
        case 5
            if norm(Ce1(ne,:))>0
                saveCell{1}=[saveCell{1}, 9];
                savedIndex{1} = [savedIndex{1}, in1(ne)];
            end
            if norm(Ce2(nw,:))>0
                saveCell{2}=[saveCell{2}, 7];
                savedIndex{2} = [savedIndex{2}, in2(nw)];
            end
            if norm(Ce3(se,:))>0
                saveCell{3}=[saveCell{3}, 3];
                savedIndex{3} = [savedIndex{3}, in3(se)];
            end
            if norm(Ce4(sw,:))>0
                saveCell{4}=[saveCell{4}, 1];
                savedIndex{4} = [savedIndex{4}, in4(sw)];
            end
            if norm(Ce1(center,:))>0
                saveCell{1}=[saveCell{1},5];
                savedIndex{1} = [savedIndex{1}, in1(center)];
            end
            if norm(Ce1(east,:))>0
                saveCell{1}=[saveCell{1},6];
                savedIndex{1} = [savedIndex{1}, in1(east)];
            end
            if norm(Ce1(north,:))>0
                saveCell{1}=[saveCell{1},8];
                savedIndex{1}=[savedIndex{1}, in1(north)];
            end
            if norm(Ce2(west,:))>0
                saveCell{2}=[saveCell{2},4];
                savedIndex{2}=[savedIndex{2}, in2(west)];
            end
            if norm(Ce2(center,:))>0
                saveCell{2}=[saveCell{2},5];
                savedIndex{2}=[savedIndex{2},in2(center)];
            end
            if norm(Ce2(north,:))>0
                saveCell{2}=[saveCell{2},8];
                savedIndex{2}=[savedIndex{2},in2(north)];
            end
            if norm(Ce3(south,:))>0
                saveCell{3}=[saveCell{3},2];
                savedIndex{3}=[savedIndex{3},in3(south)];
            end
            if norm(Ce3(center,:))>0
                saveCell{3}=[saveCell{3},5];
                savedIndex{3}=[savedIndex{3},in3(center)];
            end
            if norm(Ce3(east,:))>0
                saveCell{3}=[saveCell{3},6];
                savedIndex{3}=[savedIndex{3},in3(east)];
            end
            if norm(Ce4(west,:))>0
                saveCell{4}=[saveCell{4},4];
                savedIndex{4}=[savedIndex{4},in4(west)];
            end
            if norm(Ce4(center,:))>0
                saveCell{4}=[saveCell{4},5];
                savedIndex{4}=[savedIndex{4},in4(center)];
            end
            if norm(Ce4(south,:))>0
                saveCell{4}=[saveCell{4},2];
                savedIndex{4}=[savedIndex{4},in4(south)];
            end
            
            Ce1(ne,:) = newC1(ne,:);
            Ce2(nw,:) = newC2(nw,:);
            Ce3(se,:) = newC3(se,:);
            Ce4(sw,:) = newC4(sw,:);
            Ce1(center,:) = newC1(center,:);
            Ce1(east,:) = newC1(east,:);
            Ce1(north,:) = newC1(north,:);
            Ce2(west,:) = newC2(west,:);
            Ce2(center,:) = newC2(center,:);
            Ce2(north,:) = newC2(north,:);
            Ce3(south,:) = newC3(south,:);
            Ce3(center,:) = newC3(center,:);
            Ce3(east,:) = newC3(east,:);
            Ce4(west,:) = newC4(west,:);
            Ce4(center,:) = newC4(center,:);
            Ce4(south,:) = newC4(south,:);
          
    end
end    

%over-write the element_nod indices 
in1 = newBasisSet{1};
in2 = newBasisSet{2};
in3 = newBasisSet{3};
in4 = newBasisSet{4};


%add the saved basis functions
Ce_extra = cell(1,4);
for save_index = 1:4 
    for corner_index = 1:length(saveCell{save_index})        
        switch saveCell{save_index}(corner_index)
            case 1  %sw corner  
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(sw,:)];   
            case 2   %south
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(south,:)];
            case 3  %se corner
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(se,:)];
            case 4  %west
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(west,:)];
            case 5  %center
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(center,:)];
            case 6  %east
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(east,:)];
            case 7  %nw corner
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(nw,:)];
            case 8  %north
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(north,:)];
            case 9  %ne corner
                Ce_extra{save_index} = [Ce_extra{save_index}; Ce_save{save_index}(ne,:)];
        end
    end
end
        

%return the outputs
Ce1 = [Ce1; Ce_extra{1}];
Ce2 = [Ce2; Ce_extra{2}];
Ce3 = [Ce3; Ce_extra{3}];
Ce4 = [Ce4; Ce_extra{4}];

in1 = [in1, savedIndex{1}];
in2 = [in2, savedIndex{2}];
in3 = [in3, savedIndex{3}];
in4 = [in4, savedIndex{4}];

