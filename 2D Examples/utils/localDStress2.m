function [positive_elast]=localDStress2(e1,e2,e3,lambda,miu)
% computes the positive elastic strain energy and decomposed stress


M=sqrt((e1-e2)^2+4*e3^2);
lambda1=0.5*(e1+e2)+0.5*M;
lambda2=0.5*(e1+e2)-0.5*M;

% v1=[
%     0.5+0.5*(e1-e2)*M^-1;
%     0.5-0.5*(e1-e2)*M^-1;
%     2*e3*M^-1
%     ];
% 
% v2=[
%     0.5-0.5*(e1-e2)*M^-1;
%     0.5+0.5*(e1-e2)*M^-1;
%     -2*e3*M^-1
%     ];
% 
% d2la1de2=zeros(3,3);
% d2la1de2(1,1)=0.5*M^-1-0.5*(e1-e2)^2*M^-3;
% d2la1de2(1,2)=-d2la1de2(1,1);
% d2la1de2(1,3)=-2*(e1-e2)*e3*M^-3;
% d2la1de2(2,2)=0.5*M^-1-0.5*(e1-e2)^2*M^-3;
% d2la1de2(2,3)=2*(e1-e2)*e3*M^-3;
% d2la1de2(3,3)=2*M^-1-8*e3^2*M^-3;
% 
% d2la1de2=(d2la1de2'+d2la1de2)-diag(diag(d2la1de2));
% 
% %%%Hughes
% d2Edla1_2=lambda*(H((lambda1+lambda2))*deg+H(-(lambda1+lambda2)))+2*miu*(H(lambda1)*deg+H(-lambda1));
% d2Edla2_2=lambda*(H((lambda1+lambda2))*deg+H(-(lambda1+lambda2)))+2*miu*(H(lambda2)*deg+H(-lambda2));
% d2Edla1la2=lambda*(H(lambda1+lambda2)*deg+H(-(lambda1+lambda2)));
% dEdla1=lambda*(lambda1+lambda2)*(H((lambda1+lambda2))*deg+H(-(lambda1+lambda2)))+2*miu*lambda1*(H((lambda1))*deg+H(-(lambda1)));
% dEdla2=lambda*(lambda1+lambda2)*(H((lambda1+lambda2))*deg+H(-(lambda1+lambda2)))+2*miu*lambda2*(H((lambda2))*deg+H(-(lambda2)));
% %%%
% %%% MIehe
% % % d2Edla1_2=lambda*(deg^0.5*H(lambda1)+H(-lambda1))^2+2*miu*(H(lambda1)*deg+H(-lambda1));
% % % d2Edla2_2=lambda*(deg^0.5*H(lambda2)+H(-lambda2))^2+2*miu*(H(lambda2)*deg+H(-lambda2));
% % % d2Edla1la2=lambda*(deg^0.5*H(lambda2)+H(-lambda2))*(deg^0.5*H(lambda1)+H(-lambda1));
% % % dEdla1=lambda*((deg^0.5*H(lambda1)+H(-lambda1))*lambda1+(deg^0.5*H(lambda2)+H(-lambda2))*lambda2)*(deg^0.5*H(lambda1)+H(-lambda1))+2*miu*lambda1*(H((lambda1))*deg+H(-(lambda1)));
% % % dEdla2=lambda*((deg^0.5*H(lambda1)+H(-lambda1))*lambda1+(deg^0.5*H(lambda2)+H(-lambda2))*lambda2)*(deg^0.5*H(lambda2)+H(-lambda2))+2*miu*lambda2*(H((lambda2))*deg+H(-(lambda2)));
% %%%%
% TM=v1*d2Edla1_2*v1'+d2Edla1la2*(v1*v2'+v2*v1')+...
%     v2*d2Edla2_2*v2'+(dEdla1-dEdla2)*d2la1de2;
% 
% 
% %%%%%
% stress=(dEdla1*v1+dEdla2*v2);
%Hughes
positive_elast=(0.5*lambda*(lambda1+lambda2)^2*H(lambda1+lambda2)+miu*(lambda1^2*H(lambda1)+lambda2^2*H(lambda2)));
%Miehe
% positive_elast=(0.5*lambda*(lambda1+lambda2)^2*H(lambda1+lambda2)+miu*(lambda1^2*H(lambda1)+lambda2^2*H(lambda2)));

% % % MM=sqrt(e3^2/(e1-lambda1)^2+1);
% % % Av1=[-e3/(e1-lambda1)]*MM^-1;
% % % 
% % % MM=sqrt(e3^2/(e1-lambda2)^2+1);
% % % Av2=[-e3/(e1-lambda2)]*MM^-1;
% % % 
% % % P=[Av1',Av2'];
% % % strain=P*[lambda1*(H((lambda1))*deg+H(-(lambda1))) 0;0 lambda2*(H((lambda2))*deg+H(-(lambda2)))]*P^-1;
% % % strain=[strain(1) strain(2,2) strain(1,2)]';
end

function r=H(a)

if a>=0
    r=1;
else
    r=0;
end
end