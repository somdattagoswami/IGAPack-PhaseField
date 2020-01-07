function [Bu,Bphi,D,Hu]=strainGradHu(dgdx,d2gdx2,nument,nstress,dim,kgauss,C)

Bphi = zeros(dim,nument);
Bu = zeros(nstress,dim*nument);
% Hu = zeros((dim+1),nument);
Hu = zeros(dim*6,nument);
for inode=1:nument
    % Calculation of Bu
    Bu(1,2*inode-1)=dgdx(kgauss,1,inode);
    Bu(2,2*inode)=dgdx(kgauss,2,inode);
    Bu(3,2*inode-1)=dgdx(kgauss,2,inode);
    Bu(3,2*inode)=dgdx(kgauss,1,inode);
    % Calculation of Bphi
    Bphi(1,inode)=dgdx(kgauss,1,inode);
    Bphi(2,inode)=dgdx(kgauss,2,inode);    
end
D = C*Bu;
% for inode=1:nument
%     % Calculation of Hu
% 	Hu(1,inode) = d2gdx2(kgauss,1,inode);
%     Hu(2,inode) = d2gdx2(kgauss,3,inode);
%     Hu(3,inode) = d2gdx2(kgauss,2,inode);
% end
for inode=1:nument
    Hu(1,inode) = d2gdx2(kgauss,1,inode);
    Hu(4,inode) = d2gdx2(kgauss,3,inode);
    Hu(5,inode) = d2gdx2(kgauss,3,inode);
    Hu(6,inode) = d2gdx2(kgauss,1,inode);
    Hu(7,inode) = d2gdx2(kgauss,3,inode);
    Hu(10,inode) = d2gdx2(kgauss,2,inode);
    Hu(11,inode) = d2gdx2(kgauss,2,inode);
    Hu(12,inode) = d2gdx2(kgauss,3,inode);
end
end %end function