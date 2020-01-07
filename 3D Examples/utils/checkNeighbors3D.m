function [newBasisVert,Tjunct,RTjunct,knotUl,knotUr,knotVd,knotVu,knotWb,knotWt] = checkNeighbors3D(PHTelem,e)
% Returns information about neighboring elements, to be used by the deCasteljau3dai routine
% Collect information about new basis vertices and T junctions

newBasisVert = []; %new basis vertices
Tjunct = [];  %stores the sides that will become T junctions
RTjunct = []; %stores the sides where T junctions will be removed

% Set default values for knotUl, knotUr, knotVd, knotVu, knotWb, knotWt
knotUl = 0;
knotUr = 1;
knotVd = 0;
knotVu = 1;
knotWb = 0;
knotWt = 1;

switch length(PHTelem(e).neighbor_left)
    case 0
        newBasisVert = [newBasisVert,4];
    case 1
        Tjunct = [Tjunct,4];
    case 4
        RTjunct = [RTjunct,4];
        knotUl = PHTelem(PHTelem(e).neighbor_left(1)).vertex(1);
    otherwise
        disp(['Error on element index ', num2str(e)])
        error('Left neighbor is too refined... method not implemented')
end

switch length(PHTelem(e).neighbor_right)
    case 0
        newBasisVert = [newBasisVert, 2];
    case 1
        Tjunct = [Tjunct,2];
    case 4
        RTjunct = [RTjunct,2];
        knotUr = PHTelem(PHTelem(e).neighbor_right(1)).vertex(4);
    otherwise
        disp(['Error on element index ', num2str(e)])
        error('Right neighbor is too refined... method not implemented')
end

switch length(PHTelem(e).neighbor_front)
    case 0
        newBasisVert = [newBasisVert, 1];
    case 1
        Tjunct = [Tjunct,1];
    case 4
        RTjunct = [RTjunct,1];
        knotVd = PHTelem(PHTelem(e).neighbor_front(1)).vertex(2);
    otherwise
        disp(['Error on element index ', num2str(e)])
        error('Front neighbor is too refined... method not implemented')
end

switch length(PHTelem(e).neighbor_back)
    case 0
        newBasisVert = [newBasisVert,3];
    case 1
        Tjunct = [Tjunct,3];
    case 4
        RTjunct = [RTjunct,3];
        knotVu = PHTelem(PHTelem(e).neighbor_back(1)).vertex(5);
    otherwise
        disp(['Error on element index ', num2str(e)])
        error('Back neighbor is too refined... method not implemented')
end

switch length(PHTelem(e).neighbor_down)
    case 0
        newBasisVert = [newBasisVert,5];
    case 1
        Tjunct = [Tjunct,5];
    case 4
        RTjunct = [RTjunct,5];
        knotWb = PHTelem(PHTelem(e).neighbor_down(1)).vertex(3);
    otherwise
        disp(['Error on element index ', num2str(e)])
        error('Down neighbor is too refined... method not implemented')
end

switch length(PHTelem(e).neighbor_up)
    case 0
        newBasisVert = [newBasisVert, 6];
    case 1
        Tjunct = [Tjunct,6];
    case 4
        RTjunct = [RTjunct,6];
        knotWt = PHTelem(PHTelem(e).neighbor_up(1)).vertex(6);
    otherwise
        disp(['Error on element index ', num2str(e)])
        error('Up neighbor is too refined... method not implemented')
end

switch length(PHTelem(e).neighbor_up_left)
    case 0
        if ismember(6,newBasisVert) && ismember(4,newBasisVert)
            newBasisVert = [newBasisVert, 8];
        elseif length(PHTelem(e).neighbor_up)==4 && PHTelem(e).vertex(1)==0
            RTjunct = [RTjunct, 8];
            knotWt = PHTelem(PHTelem(e).neighbor_up(1)).vertex(6);
        elseif length(PHTelem(e).neighbor_left)==4 && PHTelem(e).vertex(6)==1
            RTjunct = [RTjunct, 8];
            knotUl = PHTelem(PHTelem(e).neighbor_left(1)).vertex(1);
        end
    case 1
        Tjunct = [Tjunct, 8];
    case 2
        if (ismember(6,newBasisVert) || ismember(6,RTjunct)) && (ismember(4,newBasisVert) || ismember(4,RTjunct))
            RTjunct = [RTjunct, 8];
            knotUl = PHTelem(PHTelem(e).neighbor_up_left(1)).vertex(1);
            knotWt = PHTelem(PHTelem(e).neighbor_up_left(1)).vertex(6);
        end
    otherwise
        disp(['Error on element index ', num2str(e)])
        error('Up_left neighbor is too refined... method not implemented')
end

switch length(PHTelem(e).neighbor_down_left)
    case 0
        if ismember(5,newBasisVert) && ismember(4,newBasisVert)
            newBasisVert = [newBasisVert, 9];
        elseif length(PHTelem(e).neighbor_down)==4 && PHTelem(e).vertex(1)==0
            RTjunct = [RTjunct, 9];
            knotWb = PHTelem(PHTelem(e).neighbor_down(1)).vertex(3);
        elseif length(PHTelem(e).neighbor_left)==4 && PHTelem(e).vertex(3)==0
            RTjunct = [RTjunct, 9];
            knotUl = PHTelem(PHTelem(e).neighbor_left(1)).vertex(1);
        end
    case 1
        Tjunct = [Tjunct, 9];
    case 2
        if (ismember(5,newBasisVert) || ismember(5,RTjunct)) && (ismember(4,newBasisVert) || ismember(4,RTjunct))
            RTjunct = [RTjunct, 9];
            knotUl = PHTelem(PHTelem(e).neighbor_down_left(1)).vertex(1);
            knotWb = PHTelem(PHTelem(e).neighbor_down_left(1)).vertex(3);
        end
    otherwise
        disp(['Error on element index ', num2str(e)])
        error('Down_left neighbor is too refined... method not implemented')
end

switch length(PHTelem(e).neighbor_up_right)
    case 0
        if ismember(6,newBasisVert) && ismember(2,newBasisVert)
            newBasisVert = [newBasisVert, 10];
        elseif length(PHTelem(e).neighbor_up)==4 && PHTelem(e).vertex(4)==1
            RTjunct = [RTjunct, 10];
            knotWt = PHTelem(PHTelem(e).neighbor_up(1)).vertex(6);
        elseif length(PHTelem(e).neighbor_right)==4 && PHTelem(e).vertex(6)==1
            RTjunct = [RTjunct, 10];
            knotUr = PHTelem(PHTelem(e).neighbor_right(1)).vertex(4);
        end
    case 1
        Tjunct = [Tjunct, 10];
    case 2
        if (ismember(6,newBasisVert) || ismember(6,RTjunct)) && (ismember(2,newBasisVert) || ismember(2,RTjunct))
            RTjunct = [RTjunct, 10];
            knotUr = PHTelem(PHTelem(e).neighbor_up_right(1)).vertex(4);
            knotWt = PHTelem(PHTelem(e).neighbor_up_right(1)).vertex(6);
        end
    otherwise
        disp(['Error on element index ', num2str(e)])
        error('Up_right neighbor is too refined... method not implemented')
end

switch length(PHTelem(e).neighbor_down_right)
    case 0
        if ismember(5,newBasisVert) && ismember(2,newBasisVert)
            newBasisVert = [newBasisVert, 11];
        elseif length(PHTelem(e).neighbor_down)==4 && PHTelem(e).vertex(4)==1
            RTjunct = [RTjunct, 11];
            knotWb = PHTelem(PHTelem(e).neighbor_down(1)).vertex(3);
        elseif length(PHTelem(e).neighbor_right)==4 && PHTelem(e).vertex(3)==0
            RTjunct = [RTjunct, 11];
            knotUr = PHTelem(PHTelem(e).neighbor_right(1)).vertex(4);
        end
    case 1
        Tjunct = [Tjunct, 11];
    case 2
        if (ismember(5,newBasisVert) || ismember(5,RTjunct)) && (ismember(2,newBasisVert) || ismember(2,RTjunct))
            RTjunct = [RTjunct, 11];
            knotUr = PHTelem(PHTelem(e).neighbor_down_right(1)).vertex(4);
            knotWb = PHTelem(PHTelem(e).neighbor_down_right(1)).vertex(3);
        end
    otherwise
        disp(['Error on element index ', num2str(e)])
        error('Down_right neighbor is too refined... method not implemented')
end

switch length(PHTelem(e).neighbor_up_front)
    case 0
        if ismember(6,newBasisVert) && ismember(1,newBasisVert)
            newBasisVert = [newBasisVert, 12];
        elseif length(PHTelem(e).neighbor_up)==4 && PHTelem(e).vertex(2)==0
            RTjunct = [RTjunct, 12];
            knotWt = PHTelem(PHTelem(e).neighbor_up(1)).vertex(6);
        elseif length(PHTelem(e).neighbor_front)==4 && PHTelem(e).vertex(6)==1
            RTjunct = [RTjunct, 12];
            knotVd = PHTelem(PHTelem(e).neighbor_front(1)).vertex(2);
        end
    case 1
        Tjunct = [Tjunct, 12];
    case 2
        if (ismember(6,newBasisVert) || ismember(6,RTjunct)) && (ismember(1,newBasisVert) || ismember(1,RTjunct))
            RTjunct = [RTjunct, 12];
            knotWt = PHTelem(PHTelem(e).neighbor_up_front(1)).vertex(6);
            knotVd = PHTelem(PHTelem(e).neighbor_up_front(1)).vertex(2);
        end
    otherwise
        disp(['Error on element index ', num2str(e)])
        error('Up_front neighbor is too refined... method not implemented')
end

switch length(PHTelem(e).neighbor_down_front)
    case 0
        if ismember(5,newBasisVert) && ismember(1,newBasisVert)
            newBasisVert = [newBasisVert, 13];
        elseif length(PHTelem(e).neighbor_down)==4 && PHTelem(e).vertex(2)==0
            RTjunct = [RTjunct, 13];
            knotWb = PHTelem(PHTelem(e).neighbor_down(1)).vertex(3);
        elseif length(PHTelem(e).neighbor_front)==4 && PHTelem(e).vertex(3)==0
            RTjunct = [RTjunct, 13];
            knotVd = PHTelem(PHTelem(e).neighbor_front(1)).vertex(2);
        end
    case 1
        Tjunct = [Tjunct, 13];
    case 2
        if (ismember(5,newBasisVert) || ismember(5,RTjunct)) && (ismember(1,newBasisVert) || ismember(1,RTjunct))
            RTjunct = [RTjunct, 13];
            knotVd = PHTelem(PHTelem(e).neighbor_down_front(1)).vertex(2);
            knotWb = PHTelem(PHTelem(e).neighbor_down_front(1)).vertex(3);
        end
    otherwise
        disp(['Error on element index ', num2str(e)])
        error('Down_front neighbor is too refined... method not implemented')
end

switch length(PHTelem(e).neighbor_up_back)
    case 0
        if ismember(6,newBasisVert) && ismember(3,newBasisVert)
            newBasisVert = [newBasisVert, 14];
        elseif length(PHTelem(e).neighbor_up)==4 && PHTelem(e).vertex(5)==1
            RTjunct = [RTjunct, 14];
            knotWt = PHTelem(PHTelem(e).neighbor_up(1)).vertex(6);
        elseif length(PHTelem(e).neighbor_back)==4 && PHTelem(e).vertex(6)==1
            RTjunct = [RTjunct, 14];
            knotVu = PHTelem(PHTelem(e).neighbor_back(1)).vertex(5);
        end
    case 1
        Tjunct = [Tjunct, 14];
    case 2
        if (ismember(6,newBasisVert) || ismember(6,RTjunct)) && (ismember(3,newBasisVert) || ismember(3,RTjunct))
            RTjunct = [RTjunct, 14];
            knotVu = PHTelem(PHTelem(e).neighbor_up_back(1)).vertex(5);
            knotWt = PHTelem(PHTelem(e).neighbor_up_back(1)).vertex(6);
        end
    otherwise
        disp(['Error on element index ', num2str(e)])
        error('Up_back neighbor is too refined... method not implemented')
end

switch length(PHTelem(e).neighbor_down_back)
    case 0
        if ismember(5,newBasisVert) && ismember(3,newBasisVert)
            newBasisVert = [newBasisVert, 15];
        elseif length(PHTelem(e).neighbor_down)==4 && PHTelem(e).vertex(5)==1
            RTjunct = [RTjunct, 15];
            knotWb = PHTelem(PHTelem(e).neighbor_down(1)).vertex(3);
        elseif length(PHTelem(e).neighbor_back)==4 && PHTelem(e).vertex(3)==0
            RTjunct = [RTjunct, 15];
            knotVu = PHTelem(PHTelem(e).neighbor_back(1)).vertex(5);
        end
    case 1
        Tjunct = [Tjunct, 15];
    case 2
        if (ismember(5,newBasisVert) || ismember(5,RTjunct)) && (ismember(3,newBasisVert) || ismember(3,RTjunct))
            RTjunct = [RTjunct, 15];
            knotVu = PHTelem(PHTelem(e).neighbor_down_back(1)).vertex(5);
            knotWb = PHTelem(PHTelem(e).neighbor_down_back(1)).vertex(3);
        end
    otherwise
        disp(['Error on element index ', num2str(e)])
        error('Down_back neighbor is too refined... method not implemented')
end

switch length(PHTelem(e).neighbor_left_front)
    case 0
        if ismember(1,newBasisVert) && ismember(4,newBasisVert)
            newBasisVert = [newBasisVert, 16];
        elseif length(PHTelem(e).neighbor_front)==4 && PHTelem(e).vertex(1)==0
            RTjunct = [RTjunct, 16];
            knotVd = PHTelem(PHTelem(e).neighbor_front(1)).vertex(2);
        elseif length(PHTelem(e).neighbor_left)==4 && PHTelem(e).vertex(2)==0
            RTjunct = [RTjunct, 16];
            knotUl = PHTelem(PHTelem(e).neighbor_left(1)).vertex(1);
        end
    case 1
        Tjunct = [Tjunct, 16];
    case 2
        if (ismember(1,newBasisVert) || ismember(1,RTjunct)) && (ismember(4,newBasisVert) || ismember(4,RTjunct))
            RTjunct = [RTjunct, 16];
            knotUl = PHTelem(PHTelem(e).neighbor_left_front(1)).vertex(1);
            knotVd = PHTelem(PHTelem(e).neighbor_left_front(1)).vertex(2);
        end
    otherwise
        disp(['Error on element index ', num2str(e)])
        error('Left_front neighbor is too refined... method not implemented')
end

switch length(PHTelem(e).neighbor_right_front)
    case 0
        if ismember(1,newBasisVert) && ismember(2,newBasisVert)
            newBasisVert = [newBasisVert, 17];
        elseif length(PHTelem(e).neighbor_front)==4 && PHTelem(e).vertex(4)==1
            RTjunct = [RTjunct, 17];
            knotVd = PHTelem(PHTelem(e).neighbor_front(1)).vertex(2);
        elseif length(PHTelem(e).neighbor_right)==4 && PHTelem(e).vertex(2)==0
            RTjunct = [RTjunct, 17];
            knotUr = PHTelem(PHTelem(e).neighbor_right(1)).vertex(4);
        end
    case 1
        Tjunct = [Tjunct, 17];
    case 2
        if (ismember(2,newBasisVert) || ismember(2,RTjunct)) && (ismember(1,newBasisVert) || ismember(1,RTjunct))
            RTjunct = [RTjunct, 17];
            knotUr = PHTelem(PHTelem(e).neighbor_right_front(1)).vertex(4);
            knotVd = PHTelem(PHTelem(e).neighbor_right_front(1)).vertex(2);
        end
    otherwise
        disp(['Error on element index ', num2str(e)])
        error('Right_front neighbor is too refined... method not implemented')
end

switch length(PHTelem(e).neighbor_left_back)
    case 0
        if ismember(4,newBasisVert) && ismember(3,newBasisVert)
            newBasisVert = [newBasisVert, 18];
        elseif length(PHTelem(e).neighbor_back)==4 && PHTelem(e).vertex(1)==0
            RTjunct = [RTjunct, 18];
            knotVu = PHTelem(PHTelem(e).neighbor_back(1)).vertex(5);
        elseif length(PHTelem(e).neighbor_left)==4 && PHTelem(e).vertex(5)==1
            RTjunct = [RTjunct, 18];
            knotUl = PHTelem(PHTelem(e).neighbor_left(1)).vertex(1);
        end
    case 1
        Tjunct = [Tjunct, 18];
    case 2
        if (ismember(4,newBasisVert) || ismember(4,RTjunct)) && (ismember(3,newBasisVert) || ismember(3,RTjunct))
            RTjunct = [RTjunct, 18];
            knotUl = PHTelem(PHTelem(e).neighbor_left_back(1)).vertex(1);
            knotVu = PHTelem(PHTelem(e).neighbor_left_back(1)).vertex(5);
        end
    otherwise
        disp(['Error on element index ', num2str(e)])
        error('Left_back neighbor is too refined... method not implemented')
end

switch length(PHTelem(e).neighbor_right_back)
    case 0
        if ismember(3,newBasisVert) && ismember(2,newBasisVert)
            newBasisVert = [newBasisVert, 19];
        elseif length(PHTelem(e).neighbor_back)==4 && PHTelem(e).vertex(4)==1
            RTjunct = [RTjunct, 19];
            knotVu = PHTelem(PHTelem(e).neighbor_back(1)).vertex(5);
        elseif length(PHTelem(e).neighbor_right)==4 && PHTelem(e).vertex(5)==1
            RTjunct = [RTjunct, 19];
            knotUr = PHTelem(PHTelem(e).neighbor_right(1)).vertex(4);
        end
    case 1
        Tjunct = [Tjunct, 19];
    case 2
        if (ismember(3,newBasisVert) || ismember(3,RTjunct)) && (ismember(2,newBasisVert) || ismember(2,RTjunct))
            RTjunct = [RTjunct, 19];
            knotUr = PHTelem(PHTelem(e).neighbor_right_back(1)).vertex(4);
            knotVu = PHTelem(PHTelem(e).neighbor_right_back(1)).vertex(5);
        end
    otherwise
        disp(['Error on element index ', num2str(e)])
        error('Right_back neighbor is too refined... method not implemented')
end