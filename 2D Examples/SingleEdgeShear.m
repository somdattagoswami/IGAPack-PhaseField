% Script for a two dimensional plate under shear

close all
clear

addpath('./utils')
addpath('./example_data')
addpath('../nurbs/inst')

Input = @() Input_shearPlate;
PF.Geometry = @(geometry)initMeshTensile(geometry);
PF.Boundary = @(PHTelem,geometry,controlPts)initialBC_shear(PHTelem,geometry);
PF.History = @(gaussCord,Fract,numberElements,geometry)history_tensile(gaussCord,...
    Fract,numberElements,geometry);
PF.Trac = @(stiffUU,tdisp,Integ,dirichlet,file_name)compTreacShear(stiffUU,tdisp,...
    Integ,dirichlet,file_name);

file_name = 'FD-shear.txt';
output = fopen(file_name,'w');
fprintf(output,'%14.6e %14.6e\n',0,0);
fclose(output);

order = input('Deciding the phase field model: \n Second-order model : 1 \n Fourth-order model : 2 \n Your choice: ','s');
degrad = input('Deciding the stress degradation function: \n 1: Quadratic stress-degradation function \n 2: Cubic stress-degradation function \n Your choice: ','s');

if order == '1'
    if degrad == '1'
        
        PF.StiffUU = @(PHTelem,sizeBasis,numberElements,dgdx,shape,Fract,Mater,...
            volume,tdisp,geometry)gStiffnessUU(PHTelem,sizeBasis,...
            numberElements,dgdx,shape,Mater,volume,tdisp,geometry);
        PF.StiffPhiPhi = @(PHTelem,sizeBasis,numberElements,dgdx,shape,Fract,...
            Mater,volume,geometry,fenerg,tdisp)gStiffnessPhiPhi(PHTelem,...
            sizeBasis,numberElements,dgdx,shape,Fract,Mater,volume,geometry,fenerg);
        
        solver_2nd(Input,PF,file_name);
    else
        
        PF.StiffUU = @(PHTelem,sizeBasis,numberElements,dgdx,shape,Fract,Mater,...
            volume,tdisp,geometry)gStiffnessUUcubic(PHTelem,sizeBasis,...
            numberElements,dgdx,shape,Fract,Mater,volume,tdisp,geometry);
        PF.StiffPhiPhi = @(PHTelem,sizeBasis,numberElements,dgdx,shape,Fract,...
            Mater,volume,geometry,fenerg,tdisp)gStiffnessPhiPhicubic(PHTelem,...
            sizeBasis,numberElements,dgdx,shape,Fract,Mater,volume,geometry,fenerg,tdisp);
        
        solver_2nd(Input,PF,file_name);
    end
else
    if degrad == '1'
        
        PF.StiffUU = @(PHTelem,sizeBasis,numberElements,dgdx,shape,Fract,Mater,...
            volume,tdisp,geometry)gStiffnessUU(PHTelem,sizeBasis,...
            numberElements,dgdx,shape,Mater,volume,tdisp,geometry);
        PF.StiffPhiPhi = @(PHTelem,sizeBasis,numberElements,dgdx,d2gdx2,shape,Fract,Mater,...
            volume,fenerg,geometry,tdisp)gStiffnessPhiPhi4th(PHTelem,sizeBasis,...
            numberElements,dgdx,d2gdx2,shape,Fract,Mater,volume,fenerg,geometry);
        
        solver_4th(Input,PF,file_name);
    else
        
        PF.StiffUU = @(PHTelem,sizeBasis,numberElements,dgdx,shape,Fract,Mater,...
            volume,tdisp,geometry)gStiffnessUUcubic(PHTelem,sizeBasis,...
            numberElements,dgdx,shape,Fract,Mater,volume,tdisp,geometry);
        PF.StiffPhiPhi = @(PHTelem,sizeBasis,numberElements,dgdx,d2gdx2,shape,Fract,Mater,...
            volume,fenerg,geometry,tdisp)gStiffnessPhiPhi4thcubic(PHTelem,sizeBasis,...
            numberElements,dgdx,d2gdx2,shape,Fract,Mater,volume,fenerg,geometry,tdisp);
        
        solver_4th(Input,PF,file_name);
    end
end