function [fenerg]=history_bimaterial(gaussCord,Fract,numberElements,geometry)

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
            
            if (gaussCord{iElem}(kgauss,2) <= 20.0)
                cenerg = Fract.cenerg1;
            else
                cenerg = Fract.cenerg2;
            end
            
            if (gaussCord{iElem}(kgauss,2)>2.0)
                dis = sqrt((gaussCord{iElem}(kgauss,1)-20.0)^2+(gaussCord{iElem}(kgauss,2)-2.0)^2);
                if dis <= Fract.constl/2
                    fenerg(iElem,kgauss) = constB*cenerg*(1-dis/(Fract.constl/2))/(2*Fract.constl);
                end
            elseif (gaussCord{iElem}(kgauss,2)<=2.0)
                dis = abs((gaussCord{iElem}(kgauss,1)-20.0));
                if dis <= Fract.constl/2
                    fenerg(iElem,kgauss) = constB*cenerg*(1-dis/(Fract.constl/2))/(2*Fract.constl);
                end
            end            
        end
    end
end
end