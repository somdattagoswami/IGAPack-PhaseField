function [fenerg,elemIndex]=history_penny(gaussCord,Fract,PHTelem,geometry,numberElements)

nGaussX = geometry.ngaussX;
nGaussY = geometry.ngaussY;
nGaussZ = geometry.ngaussZ;
constB = geometry.B;
fenerg = cell(1,numberElements); % Initilaizing the Elastic Strain Energy Matrix
elemIndex = cell(1,numberElements);
elementCounter = 0;

norm_vec = [-cosd(45); 0; cosd(45)];
radius = 0.005;
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
                        
                        mag_gauss = gaussCord{elementCounter}(:,kgauss)'*norm_vec;
                        proj_gauss = (gaussCord{elementCounter}(:,kgauss) - mag_gauss*norm_vec)';
                        
                        if (norm(proj_gauss) <= radius)
                            
                            vec3 = gaussCord{elementCounter}(:,kgauss)' - proj_gauss;
                            dis = norm(vec3);
                            
                        else
                            
                            unitvec1 = proj_gauss/norm(proj_gauss);
                            vec2 = radius*unitvec1;
                            vec3 = gaussCord{elementCounter}(:,kgauss)' - vec2;
                            dis = norm(vec3); 
                            
                        end
                        if dis <= Fract.constl/2
                            localfenerg(kgauss) = constB*Fract.cenerg*(1-dis/(Fract.constl/2))/(2*Fract.constl);
                        end
                    end
                end
            end
            fenerg{elementCounter} = localfenerg;
        end
    end
end
end
