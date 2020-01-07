function cellDataNew = blockSlice(cellData, blockIndex)
%partitions the cell Data into blocks (slices) suitable for parfor

numBlocks = length(blockIndex);
cellDataNew = cell(1,numBlocks);

for iBlock =1:numBlocks    
    cellDataNew{iBlock} = cellData(blockIndex{iBlock});    
end
