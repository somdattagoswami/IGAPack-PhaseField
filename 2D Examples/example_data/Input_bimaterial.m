function [geometry, Mater, Fract, Integ] = Input_bimaterial
% Input File for two dimensional plate made of 2 materials subjected to tension
% Use C1 cubics with adaptivity

geometry = struct;
geometry.dim = 2;
geometry.nstress = 3; % Number of stress components.
geometry.patchBoundaries = {};

% Consider cubic B-Splines
geometry.p = 3; % Degree of the polynomial in U Direction
geometry.q = 3; % Degree of the polynomial in V Direction
geometry.L = 40; % Length of the plate in the X direction
geometry.W = 40; % Width of the plate in the Y direction
geometry.numPatches = 1; % Total number of patches
geometry.numElemU = 6; % Initial mesh
geometry.numElemV = 6;
geometry.toler = 1e-5;
geometry.ngaussX = geometry.p+1;
geometry.ngaussY = geometry.q+1;
geometry.maxRefLevel = 6;
geometry.threshPhi = 0.5; % Threshold for Refinement 
geometry.B = 1e3;% Parameter for initial history function

% Material properties
Mater = struct;
Mater.E1 = 37.7*1e3; % Young's Modulus based on (N/mm2)
Mater.E2 = 377*1e3; % Young's Modulus based on (N/mm2)
Mater.nu = 0.2; % Poisson's Ratio
% Mater.C = Mater.E/(1-Mater.nu^2)*[1, Mater.nu, 0; Mater.nu, 1, 0; 0, 0, (1-Mater.nu)/2]; % Plane Stress
Mater.C1 = (Mater.E1/((1+Mater.nu)*(1-2*Mater.nu)))*[ 1-Mater.nu,Mater.nu,0;Mater.nu,1-Mater.nu,0;0,0,0.5-Mater.nu]; % Plane Strain
Mater.lamda1 = Mater.nu*Mater.E1/((1+Mater.nu)*(1-2*Mater.nu)); % Lame Constant
Mater.mu1 = Mater.E1/(2*(1+Mater.nu)); % Lame Constant
Mater.C2 = (Mater.E2/((1+Mater.nu)*(1-2*Mater.nu)))*[ 1-Mater.nu,Mater.nu,0;Mater.nu,1-Mater.nu,0;0,0,0.5-Mater.nu]; % Plane Strain
Mater.lamda2 = Mater.nu*Mater.E2/((1+Mater.nu)*(1-2*Mater.nu)); % Lame Constant
Mater.mu2 = Mater.E2/(2*(1+Mater.nu)); % Lame Constant

% Properties for Fracture
Fract = struct;
Fract.cenerg1 = 1; % Critical energy release for unstable crack (Gc)
Fract.cenerg2 = 10; % Critical energy release for unstable crack (Gc)
Fract.constl = 0.3; % L0 : Length parameter which controls the spread of the damage
Fract.s = 1e-4; % Cubic degradation parameter

% Time Integration Paramenters
Integ = struct;
Integ.nstep=20000; % Number of Displacement steps
Integ.tfacto = 0; % Total increment factor for displacement increment
Integ.dfacto1 = 1e-4; % Displacement increment per time steps upto numStepsLimit
Integ.dfacto2 = 1e-5; % Displacement increment per time steps
Integ.numStepsLimit = 400;
Integ.nprint = 5; % Printing the results after how many steps
end