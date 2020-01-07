function [sol0] = applyBoundary3D(dirichlet,tfacto,stiffness,tdisp,RHS)
%apply boundary conditions while maintaining symmetry and solve the linear system
%apply full displacement
dim = 3;
fixedPts = length(dirichlet.XYZ);
bcdof = [];
bcval = [];
if nargin<5
    RHS = zeros(size(tdisp));
end
for ivfix = 1:fixedPts
    lnode = dirichlet.XYZ(ivfix);
    for idofn =1:dim
        if (dirichlet.ValXYZ(ivfix,idofn) == 1)
            itotv =(lnode-1)*dim +idofn;                      
            dispFull = dirichlet.restrainedPts(ivfix,idofn).*tfacto;
            bcdof = [bcdof,itotv];
            bcval = [bcval,dispFull];
        end
    end
end

tic
dof_all = 1:length(tdisp);
dof_int = setdiff(dof_all,bcdof);
stiffness2 = stiffness(dof_int,dof_int);
RHS=RHS-stiffness(:,bcdof)*bcval';
RHS2 = RHS(dof_int);
toc

% sol0_int = stiffness2\RHS2;
times = size(stiffness2,1);
alpha = max(sum(abs(stiffness2),2)./diag(stiffness2))-2;
%alpha = 1e4;
L1 = ichol(stiffness2, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
[sol0_int,fl1,rr1,it1] = pcg(stiffness2,RHS2,1e-4,times,L1,L1');

fprintf('PCG exited with flag %d\n', fl1)
fprintf('Residual value: %1.15g\n', rr1)
fprintf('Number of iterations: %d\n', it1)
sol0 = zeros(length(tdisp),1);
sol0(bcdof) = bcval;
sol0(dof_int) = sol0_int;

end %endfunction