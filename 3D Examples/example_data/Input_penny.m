function [geometry, Mater, Fract, Integ] = Input_penny
% Input File for cube with a penny shaped crack in the center

% Use C1 cubics with adaptivity
geometry = struct;
geometry.dim = 3;
geometry.nstress = 6; % Number of stress components.
geometry.patchBoundaries = {};

% Consider cubic B-Splines
geometry.p = 3; % Degree of the polynomial in U Direction
geometry.q = 3; % Degree of the polynomial in V Direction
geometry.r = 3;% Degree of the polynomial in W Direction
geometry.L = 0.02; % Length of the plate in the X direction
geometry.W = 0.02; % Width of the plate in the Y direction
geometry.H = 0.02; % Height of the plate in Z Direction
geometry.numPatches = 1; % Total number of patches
geometry.numElemU = 2; % Initial mesh
geometry.numElemV = 2;
geometry.numElemW = 2; 
geometry.toler = 1e-2;
geometry.ngaussX = geometry.p+1;
geometry.ngaussY = geometry.q+1;
geometry.ngaussZ = geometry.r+1;
geometry.maxRefLevel = 4;
geometry.threshPhi = 0.5; % Threshold for Refinement 
geometry.B = 1e3;% Parameter for initial history function

% Material properties
Mater = struct;
Mater.E = 20.8*1e3; % Young's Modulus based on (N/mm2)(Miehe's Paper)
Mater.nu = 0.3; % Poisson's Ratio
Mater.C = zeros(6,6);
Mater.C(1:3,1:3)= Mater.E/(1+Mater.nu)/(1-2*Mater.nu)*[1-Mater.nu,Mater.nu,Mater.nu; Mater.nu,1-Mater.nu,Mater.nu; Mater.nu,Mater.nu,1-Mater.nu];
Mater.C(4:6,4:6)= Mater.E/(1+Mater.nu)*eye(3)/2;
Mater.lamda = Mater.nu*Mater.E/((1+Mater.nu)*(1-2*Mater.nu)); % Lame Constant
Mater.mu = Mater.E/(2*(1+Mater.nu)); % Lame Constant

% Properties for Fracture
Fract = struct;
Fract.cenerg = 0.5; % Critical energy release for unstable crack (Gc)
Fract.constl = 1/(16*50); % L0 : Length parameter which controls the spread of the damage
Fract.s = 1e-4; % Cubic degradation parameter

% Time Integration Paramenters
Integ = struct;
Integ.nstep=10000; % Number of Displacement steps
Integ.tfacto = 0; % Total increment factor for displacement increment
Integ.dfacto1 = 1e-4; % Displacement increment per time steps upto numStepsLimit
Integ.dfacto2 = 1e-5; % Displacement increment per time steps
Integ.numStepsLimit = 13;
Integ.nprint = 1; % Printing the results after how many steps
end