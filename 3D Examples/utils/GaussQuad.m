function [gp,gw]=GaussQuad(ngauss)
% Returns the Gauss points and Gauss weights

if (ngauss==2)
    gp(1) = -1/sqrt(3);
    gp(2) = 1/sqrt(3);
    gw(1) = 1;
    gw(2) = 1;
    
elseif (ngauss==3)
    gp(1) = -sqrt(3/5);
    gp(2) = 0;
    gp(3) = sqrt(3/5);
    
    gw(1) = 5/9;
    gw(2) = 8/9;
    gw(3) = 5/9;
elseif (ngauss==4)
    gp(1) = -0.86113631159405257524;
    gp(2) = -0.33998104358485626481;
    gp(3) = 0.33998104358485626481;
    gp(4) = 0.86113631159405257524;
    gw(1) = 0.34785484513745385736;
    gw(2) = 0.65214515486254614264;
    gw(3) = 0.65214515486254614264;
    gw(4) = 0.34785484513745385736;
elseif(ngauss==5)
    gp(1) = -.90617984593866399282;
    gp(2) = -.53846931010568309105;
    gp(3) = 0.0;
    gp(4) = .53846931010568309105;
    gp(5) = .90617984593866399282;
    gw(1) = .23692688505618908749;
    gw(2) = .47862867049936646808;
    gw(3) = .56888888888888888888;
    gw(4) = .47862867049936646808;
    gw(5) = .23692688505618908749;
elseif(ngauss==6)
    gp(1) = -.9324695142031520;
    gp(2) = -.6612093864662645;
    gp(3) = -.2386191860831969;
    gp(4) = .2386191860831969;
    gp(5) = .6612093864662645;
    gp(6) = .9324695142031520;
    gw(1) = .1713244923791703;
    gw(2) = .3607615730481386;
    gw(3) = .4679139345726911;
    gw(4) = .4679139345726911;
    gw(5) = .3607615730481386;
    gw(6) = .1713244923791703;
elseif(ngauss==7)
    gp(1) = -.9491079123427585;
    gp(2) = -.7415311855993944;
    gp(3) = -.4058451513773972;
    gp(4) = 0.0;
    gp(5) = .4058451513773972;
    gp(6) = .7415311855993944;
    gp(7) = .9491079123427585;
    gw(1) = .1294849661688697;
    gw(2) = .2797053914892767;
    gw(3) = .3818300505051189;
    gw(4) = .4179591836734694;
    gw(5) = .3818300505051189;
    gw(6) = .2797053914892767;
    gw(7) = .1294849661688697;
elseif(ngauss==8)
    gp(1) = -.9602898564975362;
    gp(2) = -.7966664774136267;
    gp(3) = -.5255324099163290;
    gp(4) = -.1834346424956498;
    gp(5) = .1834346424956498;
    gp(6) = .5255324099163290;
    gp(7) = .7966664774136267;
    gp(8) = .9602898564975362;
    gw(1) = .1012285362903763;
    gw(2) = .2223810344533745;
    gw(3) = .3137066458778873;
    gw(4) = .3626837833783620;
    gw(5) = .3626837833783620;
    gw(6) = .3137066458778873;
    gw(7) = .2223810344533745;
    gw(8) = .1012285362903763;
elseif(ngauss==9)
    gp(1) = -.9681602395076261;
    gp(2) = -.8360311073266358;
    gp(3) = -.6133714327005904;
    gp(4) = -.3242534234038089;
    gp(5) = 0.0;
    gp(6) = .3242534234038089;
    gp(7) = .6133714327005904;
    gp(8) = .8360311073266358;
    gp(9) = .9681602395076261;
    gw(1) = .0812743883615744;
    gw(2) = .1806481606948574;
    gw(3) = .2606106964029354;
    gw(4) = .3123470770400029;
    gw(5) = .3302393550012598;
    gw(6) = .3123470770400028;
    gw(7) = .2606106964029355;
    gw(8) = .1806481606948574;
    gw(9) = .0812743883615744;
elseif(ngauss==10)
    gp(1) = -0.97390653;
    gp(2) = -0.86506337;
    gp(3) = -0.67940957;
    gp(4) = -0.43339539;
    gp(5) = -0.14887434;
    gp(6) = 0.14887434;
    gp(7) = 0.43339539;
    gp(8) = 0.67940957;
    gp(9) = 0.86506337;
    gp(10) = 0.97390653;
    gw(1) = 0.06667134;
    gw(2) = 0.14945135;
    gw(3) = 0.21908636;
    gw(4) = 0.26926672;
    gw(5) = 0.29552422;
    gw(6) = 0.29552422;
    gw(7) = 0.26926672;
    gw(8) = 0.21908636;
    gw(9) = 0.14945135;
    gw(10) = 0.06667134;
else %(ngauss==15)
    gp(1) = -.9879925180204854;
    gp(2) = -.9372733924007059;
    gp(3) = -.8482065834104272;
    gp(4) = -.7244177313601700;
    gp(5) = -.5709721726085388;
    gp(6) = -.3941513470775634;
    gp(7) = -.2011940939974345;
    gp(8) = 0.0;
    gp(9) = .2011940939974345;
    gp(10) = .3941513470775634;
    gp(11) = .5709721726085388;
    gp(12) = .7244177313601700;
    gp(13) = .8482065834104272;
    gp(14) = .9372733924007059;
    gp(15) = .9879925180204854;
    gw(1) = .03075324199611807;
    gw(2) = .07036604748811134;
    gw(3) = .1071592204671351;
    gw(4) = .1395706779261761;
    gw(5) = .1662692058169852;
    gw(6) = .1861610000155741;
    gw(7) = .1984314853271374;
    gw(8) = .2025782419255562;
    gw(9) = .1984314853271374;
    gw(10) = .1861610000155741;
    gw(11) = .1662692058169852;
    gw(12) = .1395706779261761;
    gw(13) = .1071592204671351;
    gw(14) = .07036604748811134;
    gw(15) = .03075324199611807;
end