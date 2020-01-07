function [sol] = applyBoundary(stiffness,RHS,bcdof,bcval)

bcwt=mean(diag(stiffness)); % a measure of the average  size of an element in K
RHS=RHS-stiffness(:,bcdof)*bcval';  % modify the  force vector
RHS(bcdof) = bcwt*bcval;
stiffness(bcdof,:)=0;  % zero out the rows and  columns of the K matrix
stiffness(:,bcdof) = 0;
stiffness(bcdof,bcdof)=bcwt*speye(length(bcdof));  % put ones*bcwt on the diagonal
sol = stiffness\RHS;
end


