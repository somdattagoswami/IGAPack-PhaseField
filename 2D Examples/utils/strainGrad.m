function [Bu,Bphi,D]=strainGrad(dgdx,nument,nstress,dim,kgauss,C)

Bphi = zeros(dim,nument);
Bu = zeros(nstress,dim*nument);

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

end %end function