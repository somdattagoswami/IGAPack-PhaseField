function blockIndex = blockDecompose(numberElements,elementBlockSize)
%partitions the set 1:numberElements into disjoint subsets with up to
%elementBlockSize elements


elementCounter = 0;
numBlocks = ceil(numberElements/elementBlockSize);
blockIndex = cell(1, numBlocks);
indexCounter = 0;

while elementCounter<numberElements
    indexCounter = indexCounter + 1;
    if elementCounter+elementBlockSize <= numberElements
        blockIndex{indexCounter} = elementCounter+1:elementCounter+elementBlockSize;        
    else
        blockIndex{indexCounter} = elementCounter+1:numberElements;        
    end
    elementCounter = elementCounter + elementBlockSize;
end




