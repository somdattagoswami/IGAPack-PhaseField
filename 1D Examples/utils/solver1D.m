function convergence = solver1D(PF,Deg,Sol,Input)
% Solver for 1D problems, replaces main file

% Read model and solver parameters from the input file
[Mesh,Mater,Fract] = Input();
%Preprocessing
[PHTelem,controlPts,dirichlet,Mesh] = initMesh(Mesh);
[PHTelem,Mesh.sizeBasis] = zipConforming1D(PHTelem,Mesh.dimBasis);
[Basis] = PF.CalcBasis(PHTelem,controlPts,Mesh);

%Processing
Fract.constl = 2*Fract.constl;
[fenerg] = history1D(Basis.gaussCord,Fract,PHTelem,Mesh);
Fract.constl = Fract.constl/2;
solPhi = zeros(Mesh.sizeBasis,1);
tdisp = zeros((Mesh.dim+1)*Mesh.sizeBasis,1); % Solution Vector

normInnerStep = Inf;
miter = 0;
convergence=[];
maxInnerIter = 3;
innerIter = 0;
while (normInnerStep > Mesh.toler)
    disp('Assembling the stiffness matrix.')
    tic
    [stiffUU,Rhs,elemRef] = gStiffnessUU1DMark(PHTelem,Basis,Fract,Mater,tdisp,Mesh,Deg.Fun,Sol.RHS);
    toc
    
    disp('Imposing boundary conditions and solving.')
    tic
    [solU]=applyBoundary(stiffUU,Rhs,dirichlet.X,dirichlet.ValX);
    
    disp(['Number of degrees of freedom: ', num2str(length(solU))])
    toc
    [l2relerr] = calcErrorNorms_1D(solU,PHTelem,controlPts,Mesh.p,Sol.ExactU);
    
    if (miter==0)
        miter = miter + 1;
        innerIter = 0;
    else
        innerIter = innerIter + 1;
    end
    
    convergence(miter,1:2) = [length(solU),l2relerr];
    tic
    disp('Update and internal forces.')
    tdisp(1:Mesh.sizeBasis) = solU;
    [fenerg] = internalForces1D_1Mesh(PHTelem,Basis.dgdx,tdisp,Mater,fenerg,Mesh.nGauss);
    
    toc
    disp('Updating phase field...')
    tic
    
    [stiffPhiPhi,RHSPhi] = PF.gStiffnessPhiPhi(PHTelem,Basis,Mesh,Fract,fenerg,tdisp,Deg.Deriv);
    solPhiOld = solPhi;
    solPhi = stiffPhiPhi\RHSPhi;
%     close all
    plot1 = subplot(1,2,1);
    cla(plot1)
    plotErrorSol1DIso(PHTelem,controlPts,solPhi,Mesh.p,Sol.ExactPhi,Fract)
    miter
    innerIter
    normInnerStep = norm(solPhi-solPhiOld)/sqrt(Mesh.sizeBasis)
    toc
    tic
    tdisp(Mesh.sizeBasis+1:end) = solPhi;
    refFlag = 0;
    [solPhiPatch] = transferFieldGlob2Loc(PHTelem,Mesh.sizeBasis,solPhi);
    plot2 = subplot(1,2,2);
    cla(plot2)
    plotErrorSol1DIso(PHTelem,controlPts,solU,Mesh.p,Sol.ExactU)
    for iPatch = 1:Mesh.numPatches
        if sum(elemRef{iPatch})>0
            % Set the flag to refine and update the mesh
            refFlag = 1;
        end
    end
    
    if refFlag && (innerIter == maxInnerIter)
        for iPatch = 1:Mesh.numPatches
            if sum(elemRef{iPatch})>0
                % Refine and update the mesh
                disp(['In patch ',num2str(iPatch),' refining ',num2str(sum(elemRef{iPatch})),' elements.'])
                [PHTelem{iPatch}, controlPts{iPatch}, Mesh, solPhiPatch{iPatch}] = ...
                    refineElemProjGradedIso1D(elemRef{iPatch},PHTelem{iPatch},...
                    controlPts{iPatch},Mesh,solPhiPatch{iPatch},iPatch);
            end
        end
        
        [PHTelem,Mesh.sizeBasis] = zipConforming1D(PHTelem,Mesh.dimBasis);
        [dirichlet] = initialCrack1D(PHTelem{1});
        solPhi = transferFieldLoc2Glob(PHTelem,Mesh.sizeBasis,solPhiPatch);
        tdisp = zeros((Mesh.dim+1)*Mesh.sizeBasis,1); % Solution Vector
        tdisp(Mesh.sizeBasis+1:end) = solPhi;
        [Basis]=PF.CalcBasis(PHTelem,controlPts,Mesh);
        [fenerg] = history1D(Basis.gaussCord,Fract,PHTelem,Mesh);
        normInnerStep = Inf; % Iterate at least one more time after refinement
        miter = miter + 1;
        innerIter = 0;
    end
    
    if innerIter > maxInnerIter
        break
    end
    toc
end
