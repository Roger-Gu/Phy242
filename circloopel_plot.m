clear all
close all
format long

R=10;
M=100;

rP=0;
phiP=0;
zP=25;

[Ex,Ey,Ez]=circloopel(R,M,rP,phiP,zP)
Er=sqrt(Ex^2+Ey^2)
Ez_anal=2*pi*R*zP/(R^2+zP^2)^(3/2)
%

s=0.05;
rP=0:2:20;
phiP=(0:s:1)*2*pi;
zP=25;

[RP,PHIP,ZP]=meshgrid(rP,phiP,zP);
XP=RP.*cos(PHIP);
YP=RP.*sin(PHIP);

[Ex,Ey,Ez]=circloopel(R,M,RP,PHIP,ZP);

% scale=0.75;
% quiver3(XP,YP,ZP,Ex,Ey,Ez,scale)
quiver3(XP,YP,ZP,Ex,Ey,Ez)
view(0,90)
%
% DDD