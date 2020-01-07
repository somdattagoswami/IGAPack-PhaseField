radius = 34.29;
center = [0 0 0];
sang = 0;
eang = pi/2;
arccirc = nrbcirc(radius,center,sang,eang);
arccirc.coefs = [0 0 34.29; 0 34.29/sqrt(2) 34.29; 0 0 0; 1 1/sqrt(2) 1];
U = (arccirc.coefs(1:3,3) - arccirc.coefs(1:3,1));
U = U/norm(U);
angPlane = pi/4;
RPlane = [cos(angPlane)+U(1)^2*(1-cos(angPlane)) U(1)*U(2)*(1-cos(angPlane))-U(3)*sin(angPlane) U(1)*U(3)*(1-cos(angPlane))+U(2)*sin(angPlane)
          U(2)*U(1)*(1-cos(angPlane))+U(3)*sin(angPlane) cos(angPlane)+U(2)^2*(1-cos(angPlane)) U(2)*U(3)*(1-cos(angPlane))-U(1)*sin(angPlane)
          U(3)*U(1)*(1-cos(angPlane))-U(2)*sin(angPlane) U(3)*U(2)*(1-cos(angPlane))+U(1)*sin(angPlane) cos(angPlane)+U(3)^2*(1-cos(angPlane))];
      
curv1 = arccirc.coefs(1:3,:);
curv1(1:3,2) = curv1(1:3,2)/arccirc.coefs(4,2);
curv1 = RPlane*arccirc.coefs(1:3,:);
curv1(1:3,2) = curv1(1:3,2)*arccirc.coefs(4,2);
% NurbsTest = arccirc;
arccirc.coefs(1:3,:) = curv1;

mirrorCtrlPts(1:3,1:3) = arccirc.coefs(1:3,:);
mirrorCtrlPts(1,2) = arccirc.coefs(2,2);
mirrorCtrlPts(2,2) = arccirc.coefs(1,2);
mirrorCtrlPts(1:3,2) = mirrorCtrlPts(1:3,2)/sqrt(2);

ang=-60;
Rz = [cosd(ang) -sind(ang) 0; sind(ang) cosd(ang) 0; 0 0 1];
curv2 = Rz*mirrorCtrlPts;
curv2(:,2) = curv2(:,2)*1/sqrt(2);

arccirc2 = arccirc;
arccirc2.coefs = [curv2;1 1/sqrt(2) 1];

figure
nrbctrlplot(arccirc)
hold on
nrbctrlplot(arccirc2)