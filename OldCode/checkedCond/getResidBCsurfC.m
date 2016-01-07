function [ResBCsurfC] = getResidBCsurfC(CguessL1,CguessL2,Tsurface,hguess,tglobal);
% (b) Concentration at the surface - in the documentation \tilde{G} (capital G)
% NOTE: Can be explicitly time dependent - provide time as a parameter -
% not as an element of the state vector

% depends on CguessL(2), CguessL(1), Tguess(1) - surface temperature,
% hguess, tglobal

%ResBCsurfC = CguessL1 - 1.0;
% COnst at x=0
% option (b) constant conc
%normalize by sulfer concentration in chamber
 ResBCsurfC = CguessL1(1) - 1.0;
     
%    % option (a) constant flux
%    ResBCsurfC = sulfFlux * hguess*dx1 + Dliquid *(CguessL(2) - CguessL(1));
    


% CHECK

%ResBCsurfC = 1.0 -hguess.*(CguessL2-CguessL1) +Tsurface*CguessL1;