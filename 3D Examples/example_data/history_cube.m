function [fenerg,elemIndex]=history_cube(gaussCord,Fract,PHTelem,geometry,numberElements)

nGaussX = geometry.ngaussX;
nGaussY = geometry.ngaussX;
nGaussZ = geometry.ngaussX;
constB = geometry.B;
fenerg = cell(1,numberElements); % Initilaizing the Elastic Strain Energy Matrix
elemIndex = cell(1,numberElements);
elementCounter = 0;
for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            elementCounter = elementCounter + 1;
            localfenerg = zeros(1,nGaussX*nGaussY*nGaussZ);
            elemIndex{elementCounter}(1,1:2) = [i PHTelem{indexPatch}(i).level];
            kgauss = 0;
            for ii = 1:nGaussX
                for jj = 1:nGaussY
                    for kk = 1:nGaussZ
                        kgauss = kgauss + 1;
                        if (gaussCord{elementCounter}(1,kgauss)>0.5)
                            dis = sqrt((gaussCord{elementCounter}(1,kgauss)-0.5)^2+(gaussCord{elementCounter}(3,kgauss)-0.5)^2);
                            if dis <= Fract.constl/2
                                localfenerg(kgauss) = constB*Fract.cenerg*(1-dis/(Fract.constl/2))/(2*Fract.constl);
                            end
                        elseif (gaussCord{elementCounter}(1,kgauss)<=0.5)
                            dis = abs((gaussCord{elementCounter}(3,kgauss)-0.5));
                            if dis <= Fract.constl/2
                                 localfenerg(kgauss) = constB*Fract.cenerg*(1-dis/(Fract.constl/2))/(2*Fract.constl);
                            end
                        end
                    end
                end
            end
           fenerg{elementCounter} = localfenerg; 
        end
    end
end
end
