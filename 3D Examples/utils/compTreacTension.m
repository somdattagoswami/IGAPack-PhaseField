function compTreacTension(stiffUU,tdisp,tfacto,dirichlet,file_name)
%computes the reaction forces and prints to file

dim = 3;

%establish row-index to treac array connectivity
numReact = sum(sum(dirichlet.reactForce));
fixedPts = length(dirichlet.XYZ);
rowSet = zeros(1, numReact);
rowCounter = 0;

for ivfix = 1:fixedPts
    lnode = dirichlet.XYZ(ivfix);
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
fprintf(output,'%14.6e %14.6e\n',tfacto,sumr);
fclose(output);
end

