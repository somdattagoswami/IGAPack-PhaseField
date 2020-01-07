close all
clear

addpath('./utils')
addpath('./example_data')
addpath('../../nurbs/inst')

Input = @() Input_1D;
Sol.ExactU = @(x)exact_sol_crack(x);
Sol.RHS = @(x) sin(pi*x);

DegQuadratic.Fun = @(x) (1-x).^2;
DegQuadratic.Deriv = @(x) (2);

% 2nd order phase field model
Sol.ExactPhi = @(x,Fract)exact_sol_phi2nd(x,Fract);
PF2nd.CalcBasis = @(PHTelem,controlPts,Mesh) calcBasis1D(PHTelem,controlPts,Mesh);
PF2nd.gStiffnessPhiPhi = @(PHTelem,Basis,Mesh,Fract,fenerg,tdisp,degFunDeriv) ...
    gStiffnessPhiPhi1D(PHTelem,Basis,Mesh,Fract,fenerg,tdisp,degFunDeriv);
convergence2nd = solver1D(PF2nd,DegQuadratic,Sol,Input);

% 4th order phase field model
Sol.ExactPhi = @(x,Fract)exact_sol_phi4th(x,Fract);
PF4th.CalcBasis = @(PHTelem,controlPts,Mesh) calcBasis1D4th(PHTelem,controlPts,Mesh);
PF4th.gStiffnessPhiPhi = @(PHTelem,Basis,Mesh,Fract,fenerg,tdisp,degFunDeriv) ...
    gStiffnessPhiPhi1D4th(PHTelem,Basis,Mesh,Fract,fenerg,tdisp,degFunDeriv);
convergence4th = solver1D(PF4th,DegQuadratic,Sol,Input);

% Plot the convergence graphs
figure
loglog(convergence2nd(:,1),convergence2nd(:,2),'-o','LineWidth',3)
hold on
loglog(convergence4th(:,1),convergence4th(:,2),'-o','LineWidth',3)
legend('2nd Order','4th Order')

