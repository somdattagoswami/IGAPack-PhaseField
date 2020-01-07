function solver_4th(Input,PF,file_name)

% Reads the model and solver parameters from the input file
[geometry,Mater,Fract,Integ] = Input();

elementBlockSize = 512;
addfun = @plus;
[PHTelem,controlPts,dimBasis] = PF.Geometry(geometry);
[PHTelem,sizeBasis] = zipConforming3D(PHTelem,dimBasis,geometry);

disp('Initializing boundary conditions on the initial geometry.')
[dirichlet] = PF.Boundary(PHTelem,geometry);

disp('Precomputing shape functions and derivatives.')
[shape,dgdx,d2gdx2,volume,gaussCoord,sctrxElem,numberElements]=cartdev3D2nd(PHTelem,controlPts,geometry);

disp('History function and phase field initialization.')
Fract.constl = 2*Fract.constl;
[fenerg,elemIndex] = PF.History(gaussCoord,Fract,PHTelem,geometry,numberElements);
Fract.constl = Fract.constl/2;
clear gaussCoord

disp('Domain decomposition.')
blockIndex = blockDecompose(numberElements,elementBlockSize);
numBlocks = length(blockIndex);
shapeSlice = blockSlice(shape,blockIndex);
clear shape
sctrxElemSlice = blockSlice(sctrxElem,blockIndex);
clear sctrxElem
dgdxSlice = blockSlice(dgdx,blockIndex);
clear dgdx
d2gdx2Slice = blockSlice(d2gdx2,blockIndex);
clear d2gdx2
volumeSlice = blockSlice(volume,blockIndex);
clear  volume
fenergSlice = blockSlice(fenerg,blockIndex);
clear fenerg
elemIndexSlice = blockSlice(elemIndex,blockIndex);
clear elemIndex

tdisp = zeros(4*sizeBasis,1);
solPhi = zeros(sizeBasis,1);

for istep = 1:Integ.nstep+1
    
    istep
    if (istep < Integ.numStepsLimit)
        Integ.tfacto = Integ.tfacto + Integ.dfacto1;
    else
        Integ.tfacto = Integ.tfacto + Integ.dfacto2;
    end
    
    % Begin inner iteration
    normInnerStep = Inf;
    miter = 0;
    while (normInnerStep > geometry.toler)
        refFlag = 0;
        disp('Assembling the stiffness matrix.')
        tic
        stiffUU = sparse(geometry.dim*sizeBasis,geometry.dim*sizeBasis);
        gElemIndex = cell(1,numBlocks);
        parfor iBlock = 1:numBlocks
            lengthBlock = length(blockIndex{iBlock});
            [stiffUUCell,elemRef] = PF.StiffUU(lengthBlock,sctrxElemSlice{iBlock},sizeBasis,dgdxSlice{iBlock},...
                shapeSlice{iBlock},volumeSlice{iBlock},elemIndexSlice{iBlock},Fract,Mater,geometry,solPhi);
            stiffUU = addfun(stiffUU, stiffUUCell);
            gElemIndex{iBlock} = elemRef;
        end
        clear stiffUUCell elemRef
        toc
        elemRef = zeros(1,length(PHTelem{1}));
        if sum([gElemIndex{:}])
            refFlag = 1;
            for iBlock = 1:numBlocks
                [rIndex,cIndex] = find(gElemIndex{iBlock}>0);
                elemNum = gElemIndex{iBlock}(1,cIndex);
                elemRef(1,elemNum) = 1;
            end
        end
        clear gElemIndex rIndex cIndex elemNum
        disp('Imposing boundary conditions and solving.')
        solU = applyBoundary3D(dirichlet,Integ.tfacto,stiffUU,tdisp(1:3*sizeBasis));
        
        tdisp(1:3*sizeBasis) = solU;
        clear solU
        stiffPhiPhi = sparse(sizeBasis,sizeBasis);
        RHSPhi = zeros(sizeBasis,1);
        
        disp('Update the internal forces.')
        tic
        for iBlock = 1:numBlocks
            lengthBlock = length(blockIndex{iBlock});
            [stiffPhiPhiCell,RHSPhiCell,fenergSlice{iBlock}] = PF.StiffPhiPhi(lengthBlock,...
                sctrxElemSlice{iBlock},sizeBasis,dgdxSlice{iBlock},d2gdx2Slice{iBlock},...
                shapeSlice{iBlock},volumeSlice{iBlock},Fract,Mater,geometry,fenergSlice{iBlock},...
                tdisp(1:3*sizeBasis),solPhi);
            
            stiffPhiPhi = addfun(stiffPhiPhi, stiffPhiPhiCell);
            RHSPhi = addfun(RHSPhi, RHSPhiCell);
        end
        toc
        
        solPhiOld = solPhi;
        solPhi = stiffPhiPhi\RHSPhi;
        clear stiffPhiPhi RHSPhi
        miter = miter + 1
        normInnerStep = norm(solPhi-solPhiOld)
        tdisp(3*sizeBasis+1:end) = solPhi;
        clear solPhiOld
        
        [solPhiPatch] = transferFieldGlob2Loc(PHTelem,dimBasis,solPhi);
        if refFlag
            tic
            clear stiffUU
            [PHTelem{1},controlPts{1},dimBasis,solPhiPatch{1},numberElements] = ...
                refineElemProjGradedIso3D(elemRef,PHTelem{1},controlPts{1},geometry,...
                dimBasis,solPhiPatch{1},numberElements);
            [PHTelem,sizeBasis] = zipConforming3D(PHTelem,dimBasis,geometry);
            [dirichlet] = PF.Boundary(PHTelem,geometry);
            [shape,dgdx,d2gdx2,volume,gaussCoord,sctrxElem]=cartdev3D2nd(PHTelem,controlPts,geometry);
            [fenerg,elemIndex] = PF.History(gaussCoord,Fract,PHTelem,geometry,numberElements);
            clear gaussCoord
            
            disp('Domain decomposition.')
            tic
            blockIndex = blockDecompose(numberElements,elementBlockSize);
            numBlocks = length(blockIndex);
            shapeSlice = blockSlice(shape,blockIndex);
            clear shape
            sctrxElemSlice = blockSlice(sctrxElem,blockIndex);
            clear sctrxElem
            dgdxSlice = blockSlice(dgdx,blockIndex);
            clear dgdx
            d2gdx2Slice = blockSlice(d2gdx2,blockIndex);
            clear d2gdx2
            volumeSlice = blockSlice(volume,blockIndex);
            clear  volume
            fenergSlice = blockSlice(fenerg,blockIndex);
            clear fenerg
            elemIndexSlice = blockSlice(elemIndex,blockIndex);
            clear elemIndex
            
            solPhi = transferFieldLoc2Glob(PHTelem,sizeBasis,solPhiPatch);
            tdisp = zeros(4*sizeBasis,1); % Solution Vector
            tdisp(3*sizeBasis+1:end) = solPhi;
            normInnerStep = Inf;
            toc
        end
    end
    
    disp('Print data for force-tdisp curves')
    tic
    geometry.numberElements = numberElements;
    PF.Trac(stiffUU,tdisp(1:3*sizeBasis),Integ.tfacto,dirichlet,file_name)
    if(mod(istep,Integ.nprint) == 0)% Print results
        vtuFile = ['time_',num2str(istep),'.vtu'];
        plot3D(PHTelem,tdisp,sizeBasis,geometry,controlPts,vtuFile)
        fprintf('Done step: %5d\n',istep);
        
    end %if
    toc
end %istep
end