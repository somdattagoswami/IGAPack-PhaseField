function [fenerg]=history1D(gaussCord,Fract,PHTelem,Mesh)

numPatches = length(PHTelem);
fenerg = cell(1,numPatches);

for indexPatch = 1:length(PHTelem)
    fenerg(indexPatch) = {cell(1,length(PHTelem{indexPatch}))};
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            fenerg{indexPatch}{i} = zeros(1,Mesh.nGauss);
            kgauss = 0;
            for ii=1:Mesh.nGauss
                kgauss = kgauss + 1;
                dis = abs((gaussCord{indexPatch}{i}(kgauss,1)-0));
                if dis <= Fract.constl/2
                    fenerg{indexPatch}{i}(kgauss) = Mesh.B*Fract.cenerg*(1-dis/(Fract.constl/2))/(2*Fract.constl);
                end
            end
        end
    end
end
end