function [C, nb] = bezierExtraction(knot,p)
% Bezier extraction
% Based on Algroithm 1, from Borden - Isogeometric finite element data
% structures based on Bezier extraction

m  = length(knot)-p-1;
a  = p+1;
b  = a+1;
nb = 1;
C(:,:,1) = eye(p+1);

while b <= m
    C(:,:,nb+1) = eye(p+1);
    i=b;
    while b <= m && knot(b+1) == knot(b)
        b=b+1;
    end
    
    multiplicity = b-i+1;
    if multiplicity < p
        numerator=knot(b)-knot(a);
        for j=p:-1:multiplicity+1
            alphas(j-multiplicity)=numerator/(knot(a+j)-knot(a));
        end
        r=p-multiplicity;
        for j=1:r
            save = r-j+1;
            s = multiplicity + j;
            for k=p+1:-1:s+1
                alpha=alphas(k-s);
                C(:,k,nb)=alpha*C(:,k,nb)+(1-alpha)*C(:,k-1,nb);
            end
            if b <= m
                C(save:save+j,save,nb+1)=C(p-j+1:p+1,p+1,nb);
            end
        end
        nb=nb+1;
        if b <= m
            a=b;
            b=b+1;
        end
    elseif multiplicity==p
        if b <= m
            nb=nb+1; a=b; b=b+1;
        end
    end
end

