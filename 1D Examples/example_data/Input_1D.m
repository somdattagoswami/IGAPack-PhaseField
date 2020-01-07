function [Mesh, Mater, Fract] = Input_1D

% Input File for 1D 
% Use C1 cubics with adaptivity
Mesh = struct;
Mesh.dim = 1;
Mesh.patchBoundaries = {};

% Consider cubic B-Splines
Mesh.p = 3; % Degree of the polynomial in U Direction
Mesh.L = 1; % Half the length of the beam [-L, L]
Mesh.numPatches = 1; % Total number of patches
Mesh.numElem = 12; % Initial mesh
Mesh.toler = 1e-12;
Mesh.nGauss = Mesh.p+1;
Mesh.maxRefLevel = 9;
Mesh.threshPhi = 0.2; % Threshold for Refinement 
Mesh.B = 1e3;% Parameter for initial history function

% Material properties
Mater = struct;
Mater.E = 1; % Young's Modulus based on (N/mm2)
Mater.nu = 0.3; % Poisson's Ratio

% Properties for Fracture
Fract = struct;
Fract.constk = 0; % For numerical convenience to get well conditioned system of equation
Fract.cenerg = 2.7; % Critical energy release for unstable crack (Gc)
Fract.constl = 0.0125/2; % L0 : Length parameter which controls the spread of the damage
end
