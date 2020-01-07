function [B,dB]=bernstein_basis(xi,p)
%Algorithm A1.3 in Piegl & Tiller
%xi is a scalar or a row vector
B = zeros(length(xi),p+1);
B(:,1) = 1;
u1 = 1-xi';
u2 = 1+xi';
for j=1:p
    saved = 0;
    for k=0:j-1
        temp = B(:,k+1);
        B(:,k+1) = saved + u1.*temp;
        saved = u2.*temp;
    end
    B(:,j+1) = saved;
end
B = B./(2^p);

%calculate the p-1 Bernstein polynomial for derivative calculations
dB = zeros(length(xi),p);
dB(:,1) = 1;
for j=1:p-1
    saved = 0;
    for k=0:j-1
        temp = dB(:,k+1);
        dB(:,k+1) = saved + u1.*temp;
        saved = u2.*temp;
    end
    dB(:,j+1) = saved;
end
dB = dB./(2^p);
dB = [zeros(length(xi),1), dB, zeros(length(xi),1)];
dB = (dB(:,1:end-1)-dB(:,2:end))*p;

