function [fenerg]=history_tensile(gaussCord,Fract,numberElements,geometry)

nGaussX = geometry.ngaussX;
nGaussY = geometry.ngaussX;
constB = geometry.B;
mgauss = nGaussX*nGaussY;
fenerg = zeros(numberElements,mgauss);

for iElem = 1:numberElements
    kgauss = 0;
    for ii=1:nGaussX
        for jj=1:nGaussY
            kgauss = kgauss + 1;
            if (gaussCord{iElem}(kgauss,1)>0.5)
                dis = sqrt((gaussCord{iElem}(kgauss,1)-0.5)^2+(gaussCord{iElem}(kgauss,2)-0.5)^2);
                if dis <= Fract.constl/2
                    fenerg(iElem,kgauss) = constB*Fract.cenerg*(1-dis/(Fract.constl/2))/(2*Fract.constl);
                end
            elseif (gaussCord{iElem}(kgauss,1)<=0.5)
                dis = abs((gaussCord{iElem}(kgauss,2)-0.5));
                if dis <= Fract.constl/2
                    fenerg(iElem,kgauss) = constB*Fract.cenerg*(1-dis/(Fract.constl/2))/(2*Fract.constl);
                end
            end            
        end
    end
end
end
