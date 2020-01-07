function [geometry, Mater, Fract, Integ] = Input_cencrack
% Input File for two dimensional plate with an initial crack at the center
% subjected to tensile loading

% Use C1 cubics with adaptivity
geometry = struct;
geometry.dim = 2;
geometry.nstress = 3; % Number of stress components.
geometry.patchBoundaries = {};

% Consider cubic B-Splines
geometry.p = 3; % Degree of the polynomial in U Direction
geometry.q = 3; % Degree of the polynomial in V Direction
geometry.L = 20; % Length of the plate in the X direction
geometry.W = 200; % Width of the plate in the Y direction
geometry.numPatches = 1; % Total number of patches
geometry.numElemU = 5; % Initial mesh
geometry.numElemV = 50;
geometry.toler = 0.5*1e-4;
geometry.ngaussX = geometry.p+1;
geometry.ngaussY = geometry.q+1;
geometry.maxRefLevel = 5;
geometry.threshPhi = 0.5; % Threshold for Refinement 
geometry.B = 1e3;% Parameter for initial history function

% Material properties
Mater = struct;
Mater.E = 70*1e3; % Young's Modulus based on (N/mm2)(Miehe's Paper)
Mater.nu = 0.22; % Poisson's Ratio
Mater.C = (Mater.E/((1+Mater.nu)*(1-2*Mater.nu)))*[ 1-Mater.nu,Mater.nu,0;Mater.nu,1-Mater.nu,0;0,0,0.5-Mater.nu]; % Plane Strain
Mater.lamda = Mater.nu*Mater.E/((1+Mater.nu)*(1-2*Mater.nu)); % Lame Constant
Mater.mu = Mater.E/(2*(1+Mater.nu)); % Lame Constant

% Properties for Fracture
Fract = struct;
Fract.cenerg = 0.007; % Critical energy release for unstable crack (Gc)
Fract.constl = 0.25; % L0 : Length parameter which controls the spread of the damage
Fract.s = 1e-4; % Cubic degradation parameter

% Time Integration Paramenters
Integ = struct;
Integ.nstep=10000; % Number of Displacement steps
Integ.tfacto = 0; % Total increment factor for displacement increment
Integ.dfacto1 = 2.5e-4; % Displacement increment per time steps upto numStepsLimit
Integ.dfacto2 = 2.5e-5; % Displacement increment per time steps
Integ.numStepsLimit = 35;
Integ.nprint = 5; % Printing the results after how many steps
end