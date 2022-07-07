function [Ex,Ey,Ez]=circloopel(R,M,rP,phiP,zP)
% ,phiQ_vec
% Ke=\lambda/(4*pi*epsilon0)=1;

% R=10;

% M=100;
Delta_phi=2*pi/M;
Delta_ell=R*Delta_phi;

% rP=0;
% phiP=0;
% zP=0.1;

m=0;
phiQ=0;
DeltaEx=0;
Ex=0;
DeltaEy=0;
Ey=0;
DeltaEz=0;
Ez=0;

for m=1:M

    phiQ=-Delta_phi+m*Delta_phi;
    % phiQ_vec(m)=phiQ;

    DeltaEx=Delta_ell*(rP.*cos(phiP)-R*cos(phiQ))./(rP.^2+R^2-2*R*rP.*cos(phiQ-phiP)+zP.^2).^(3/2);
    Ex=Ex+DeltaEx;

    DeltaEy=Delta_ell*(rP.*sin(phiP)-R*sin(phiQ))./(rP.^2+R^2-2*R*rP.*cos(phiQ-phiP)+zP.^2).^(3/2);
    Ey=Ey+DeltaEy;

    DeltaEz=Delta_ell*zP./(rP.^2+R^2-2*R*rP.*cos(phiQ-phiP)+zP.^2).^(3/2);
    Ez=Ez+DeltaEz;

end

end
% DDD