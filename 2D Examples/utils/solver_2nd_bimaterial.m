function solver_2nd_bimaterial(Input,PF,file_name)

% Reads the model and solver parameters from the input file
[geometry,Mater,Fract,Integ] = Input();
[PHTelem,controlPts,dimBasis] = PF.Geometry(geometry);
[PHTelem,sizeBasis] = zipConforming(PHTelem,dimBasis,geometry);

scrsz = get(groot, 'ScreenSize');
hFig = figure('Position',[1 scrsz(4)/6 3*scrsz(3)/5 3*scrsz(4)/4]);
plot1 = subplot(2,2,[1,2]);
cla(plot1)
plotMesh2D(PHTelem,controlPts,geometry)
axis equal
title('Intial Mesh');

disp('Initializing boundary conditions on the initial geometry.')
[dirichlet] = PF.Boundary(PHTelem,geometry);

disp('Precomputing shape functions and derivatives.')
[shape,dgdx,volume,gaussCord,numberElements]=cartdev(PHTelem,controlPts,geometry);

disp('History function and phase field initialization.')
Fract.constl = 4*Fract.constl;
[fenerg] = PF.History(gaussCord,Fract,numberElements,geometry);
Fract.constl = Fract.constl/4;

tdisp = zeros(3*sizeBasis,1);
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
        [stiffUU,elemRef] = PF.StiffUU(PHTelem,sizeBasis,numberElements,dgdx,...
            shape,Fract,gaussCord,Mater,volume,tdisp,geometry);
        toc
        
        disp('Imposing boundary conditions and solving.')
        tic
        solU = applyBoundary2D(dirichlet,Integ.tfacto,stiffUU,tdisp(1:2*sizeBasis));
        tdisp(1:2*sizeBasis) = solU;
        clear solU
        toc
        
        disp('Update the internal forces.')
        fenerg = internalForces_bimaterial(PHTelem,dgdx,tdisp,geometry,Mater,fenerg,gaussCord);
        
        disp('Updating phase field.')
        tic
        [stiffPhiPhi,RHSPhi] = PF.StiffPhiPhi(PHTelem,sizeBasis,numberElements,...
            dgdx,shape,Fract,Mater,volume,geometry,fenerg,gaussCord,tdisp);
        
        solPhiOld = solPhi;
        solPhi = stiffPhiPhi\RHSPhi;
        normInnerStep = norm(stiffPhiPhi*solPhiOld-RHSPhi)/norm(RHSPhi)
        miter = miter + 1
        tdisp(2*sizeBasis+1:end) = solPhi;
        clear stiffPhiPhi RHSPhi solPhiOld
        [solPhiPatch] = transferFieldGlob2Loc(PHTelem,dimBasis,solPhi);
        toc
        
        for iPatch = 1:geometry.numPatches
            if sum(elemRef{iPatch})>0
                % Refine and update the mesh
                refFlag = 1;
                disp(['In patch ',num2str(iPatch),' refining ',num2str(sum(elemRef{iPatch})),' elements.'])
                [PHTelem{iPatch},controlPts{iPatch},dimBasis(iPatch),solPhiPatch{iPatch}, ...
                    numberElements] = refineElemProjGradedIso(elemRef{iPatch},PHTelem{iPatch}, ...
                    controlPts{iPatch},geometry,dimBasis(iPatch),solPhiPatch{iPatch},numberElements);
            end
        end
        
        if refFlag
            tic
            clear stiffUU
            [PHTelem,sizeBasis] = zipConforming(PHTelem,dimBasis,geometry);
            plot1 = subplot(2,2,[1,2]);
            cla(plot1)
            plotMesh2D(PHTelem,controlPts,geometry)
            axis equal
            title(['Modified Mesh for Loadstep', num2str(istep) ,' and Iteration ', num2str(miter)]);
            [dirichlet] = PF.Boundary(PHTelem,geometry);
            [shape,dgdx,volume,gaussCord,numberElements] = cartdev(PHTelem,controlPts,geometry);
            [fenerg] = PF.History(gaussCord,Fract,numberElements,geometry);
            
            solPhi = transferFieldLoc2Glob(PHTelem,sizeBasis,solPhiPatch);
            tdisp = zeros(3*sizeBasis,1); % Solution Vector
            tdisp(2*sizeBasis+1:end) = solPhi;
            normInnerStep = Inf;
            toc
        end
    end
    
    disp('Print data for force-tdisp curves.')
    tic
    PF.Trac(stiffUU,tdisp(1:2*sizeBasis),Integ.tfacto,dirichlet,file_name);
    if(mod(istep,Integ.nprint) == 0)% Print results
        fprintf('Done step: %5d\n',istep);
        plotDispPhase2D_bimaterial(PHTelem,tdisp,sizeBasis,numberElements,geometry,controlPts,Mater,gaussCord)
        plot1 = subplot(2,2,[1,2]);
        title(['Mesh for Loadstep ',num2str(istep),' and Iteration ',num2str(miter)]);
        saveas(hFig, ['Loadstep', num2str(istep),'.png'])
    end %if
    toc    
end %istep
end