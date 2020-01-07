function compTreacTension(stiffUU,tdisp,tfacto,dirichlet,file_name)
% Computes the reaction forces and prints to file

dim = 2;
% Establish row-index to treac array connectivity
numReact = sum(sum(dirichlet.reactForce));
fixedPts = length(dirichlet.XY);
rowSet = zeros(1, numReact);
rowCounter = 0;

for ivfix = 1:fixedPts
    lnode = dirichlet.XY(ivfix);
    for idofn =1:dim
        if(dirichlet.reactForce(ivfix,idofn) == 1)
            rowCounter = rowCounter + 1;
            itotv =(lnode-1)*dim +idofn;
            rowSet(rowCounter) = itotv;
            %treac(ivfix,idofn) =treac(ivfix,idofn)- stiffness(itotv,:)*tdisp;            
        end
    end
end

treacTemp = -stiffUU(rowSet, :)*tdisp;
sumr = sum(treacTemp);

%open file for writing (append)
output = fopen(file_name,'a');
disp(['Displacement: ', num2str(tfacto), ' Reaction force: ', num2str(sumr)])
fprintf(output,'%14.6e %14.6e\n',tfacto,-sumr);
fclose(output);
end