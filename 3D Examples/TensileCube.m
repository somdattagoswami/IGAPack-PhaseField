% Script for Cube under tension

close all
clear

addpath('./utils')
addpath('./example_data')
addpath('../nurbs/inst')

Input = @() Input_cube;
PF.Geometry = @(geometry)initMeshCube(geometry);
PF.Boundary = @(PHTelem,geometry,controlPts)initialBC_cube(PHTelem,geometry);
PF.History = @(gaussCoord,Fract,PHTelem,geometry,numberElements)history_cube(gaussCoord,Fract,PHTelem,geometry,numberElements);
PF.Trac = @(stiffUU,tdisp,Integ,dirichlet,file_name)compTreacTension(stiffUU,tdisp,Integ,dirichlet,file_name);

file_name = 'FD-cube.txt';
output = fopen(file_name,'w');
fprintf(output,'%14.6e %14.6e\n',0,0);
fclose(output);

order = input('Deciding the phase field model: \n Second-order model : 1 \n Fourth-order model : 2 \n Your choice: ','s');
degrad = input('Deciding the stress degradation function: \n 1: Quadratic stress-degradation function \n 2: Cubic stress-degradation function \n Your choice: ','s');

if order == '1'
    if degrad == '1'
        
        PF.StiffUU = @(lengthBlock,sctrxElemSlice,sizeBasis,dgdxSlice,...
            shapeSlice,volumeSlice,elemIndexSlice,Fract,Mater,geometry,...
            solPhi)gStiffnessUU3D(lengthBlock,sctrxElemSlice,sizeBasis,dgdxSlice,...
            shapeSlice,volumeSlice,elemIndexSlice,Mater,geometry,solPhi);
        
        PF.StiffPhiPhi = @(lengthBlock,sctrxElemSlice,sizeBasis,dgdxSlice,...
            shapeSlice,volumeSlice,Fract,Mater,geometry,fenergSlice,tdisp,...
            solPhi)gStiffnessPhiPhi3D(lengthBlock,sctrxElemSlice,sizeBasis,dgdxSlice,...
            shapeSlice,volumeSlice,Fract,Mater,geometry,fenergSlice,tdisp);
        
        solver_2nd(Input,PF,file_name);
    else
        
        PF.StiffUU = @(lengthBlock,sctrxElemSlice,sizeBasis,dgdxSlice,...
            shapeSlice,volumeSlice,elemIndexSlice,Fract,Mater,geometry,...
            solPhi)gStiffnessUU3Dcubic(lengthBlock,sctrxElemSlice,sizeBasis,dgdxSlice,...
            shapeSlice,volumeSlice,elemIndexSlice,Fract,Mater,geometry,solPhi);
        
        PF.StiffPhiPhi = @(lengthBlock,sctrxElemSlice,sizeBasis,dgdxSlice,...
            shapeSlice,volumeSlice,Fract,Mater,geometry,fenergSlice,tdisp,...
            solPhi)gStiffnessPhiPhi3Dcubic(lengthBlock,sctrxElemSlice,sizeBasis,dgdxSlice,...
            shapeSlice,volumeSlice,Fract,Mater,geometry,fenergSlice,tdisp,solPhi);
        
        solver_2nd(Input,PF,file_name);
    end
else
    if degrad == '1'
        
        PF.StiffUU = @(lengthBlock,sctrxElemSlice,sizeBasis,dgdxSlice,...
            shapeSlice,volumeSlice,elemIndexSlice,Fract,Mater,geometry,...
            solPhi)gStiffnessUU3D(lengthBlock,sctrxElemSlice,sizeBasis,dgdxSlice,...
            shapeSlice,volumeSlice,elemIndexSlice,Mater,geometry,solPhi);
        
        PF.StiffPhiPhi = @(lengthBlock,sctrxElemSlice,sizeBasis,dgdxSlice,d2gdx2Slice,...
            shapeSlice,volumeSlice,Fract,Mater,geometry,fenergSlice,tdisp,solPhi)...
            gStiffnessPhiPhi3D4th(lengthBlock,sctrxElemSlice,sizeBasis,dgdxSlice,...
            d2gdx2Slice,shapeSlice,volumeSlice,Fract,Mater,geometry,fenergSlice,tdisp);
        
        solver_4th(Input,PF,file_name);
    else
        
        PF.StiffUU = @(lengthBlock,sctrxElemSlice,sizeBasis,dgdxSlice,...
            shapeSlice,volumeSlice,elemIndexSlice,Fract,Mater,geometry,...
            solPhi)gStiffnessUU3Dcubic(lengthBlock,sctrxElemSlice,sizeBasis,dgdxSlice,...
            shapeSlice,volumeSlice,elemIndexSlice,Fract,Mater,geometry,solPhi);
        
        PF.StiffPhiPhi = @(lengthBlock,sctrxElemSlice,sizeBasis,dgdxSlice,d2gdx2Slice,...
            shapeSlice,volumeSlice,Fract,Mater,geometry,fenergSlice,tdisp,solPhi)...
            gStiffnessPhiPhi3D4thcubic(lengthBlock,sctrxElemSlice,sizeBasis,dgdxSlice,...
            d2gdx2Slice,shapeSlice,volumeSlice,Fract,Mater,geometry,fenergSlice,tdisp,solPhi);
        
        solver_4th(Input,PF,file_name);
    end
end