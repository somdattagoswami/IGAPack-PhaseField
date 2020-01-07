%  Script for a two dimensional plate made of 2 materials subjected to tension

close all
clear

addpath('./utils')
addpath('./example_data')
addpath('../nurbs/inst')

Input = @() Input_bimaterial;
PF.Geometry = @(geometry)initMeshTensile(geometry);
PF.Boundary = @(PHTelem,geometry)initialBC_bimaterial(PHTelem,geometry);
PF.History = @(gaussCord,Fract,numberElements,geometry)history_bimaterial(gaussCord,...
    Fract,numberElements,geometry);
PF.Trac = @(stiffUU,tdisp,Integ,dirichlet,file_name)compTreacTension(stiffUU,tdisp,...
    Integ,dirichlet,file_name);

file_name = 'FD-bimaterial.txt';
output = fopen(file_name,'w');
fprintf(output,'%14.6e %14.6e\n',0,0);
fclose(output);

order = input('Deciding the phase field model: \n Second-order model : 1 \n Fourth-order model : 2 \n Your choice: ','s');
degrad = input('Deciding the stress degradation function: \n 1: Quadratic stress-degradation function \n 2: Cubic stress-degradation function \n Your choice: ','s');

if order == '1'
    if degrad == '1'
        
        PF.StiffUU = @(PHTelem,sizeBasis,numberElements,dgdx,shape,Fract,gaussCord,Mater,...
            volume,tdisp,geometry)gStiffnessUU_bimaterial(PHTelem,sizeBasis,...
            numberElements,dgdx,shape,gaussCord,Mater,volume,tdisp,geometry);
        
        PF.StiffPhiPhi = @(PHTelem,sizeBasis,numberElements,dgdx,shape,Fract,...
            Mater,volume,geometry,fenerg,gaussCord,tdisp)gStiffnessPhiPhi_bimaterial(PHTelem,...
            sizeBasis,numberElements,dgdx,shape,Fract,Mater,volume,geometry,fenerg,gaussCord);
        
        solver_2nd_bimaterial(Input,PF,file_name);
    else
        
        PF.StiffUU = @(PHTelem,sizeBasis,numberElements,dgdx,shape,Fract,gaussCord,Mater,...
            volume,tdisp,geometry)gStiffnessUUcubic_bimaterial(PHTelem,sizeBasis,...
            numberElements,dgdx,shape,gaussCord,Fract,Mater,volume,tdisp,geometry);
        
        PF.StiffPhiPhi = @(PHTelem,sizeBasis,numberElements,dgdx,shape,Fract,...
            Mater,volume,geometry,fenerg,gaussCord,tdisp)gStiffnessPhiPhicubic_bimaterial(PHTelem,...
            sizeBasis,numberElements,dgdx,shape,Fract,Mater,volume,geometry,fenerg,gaussCord,tdisp);
        
        solver_2nd_bimaterial(Input,PF,file_name);
    end
else
    if degrad == '1'
        
        PF.StiffUU = @(PHTelem,sizeBasis,numberElements,dgdx,shape,Fract,gaussCord,Mater,...
            volume,tdisp,geometry)gStiffnessUU_bimaterial(PHTelem,sizeBasis,...
            numberElements,dgdx,shape,gaussCord,Mater,volume,tdisp,geometry);
        
        PF.StiffPhiPhi = @(PHTelem,dimBasis,numberElements,dgdx,d2gdx2,shape,...
            Fract,Mater,volume,fenerg,geometry,gaussCord,tdisp)gStiffnessPhiPhi4th_bimaterial(PHTelem,...
            dimBasis,numberElements,dgdx,d2gdx2,shape,Fract,Mater,volume,fenerg,geometry,gaussCord);
        
        solver_4th_bimaterial(Input,PF,file_name);
    else
        
        PF.StiffUU = @(PHTelem,sizeBasis,numberElements,dgdx,shape,Fract,gaussCord,Mater,...
            volume,tdisp,geometry)gStiffnessUUcubic_bimaterial(PHTelem,sizeBasis,...
            numberElements,dgdx,shape,gaussCord,Fract,Mater,volume,tdisp,geometry);
        
        PF.StiffPhiPhi = @(PHTelem,dimBasis,numberElements,dgdx,d2gdx2,shape,...
            Fract,Mater,volume,fenerg,geometry,gaussCord,tdisp)gStiffnessPhiPhi4thcubic_bimaterial(PHTelem,dimBasis,...
            numberElements,dgdx,d2gdx2,shape,Fract,Mater,volume,fenerg,geometry,gaussCord,tdisp);
        
        solver_4th_bimaterial(Input,PF,file_name);
    end
end