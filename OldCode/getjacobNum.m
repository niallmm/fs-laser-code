function [jacob] = getjacobNum(uold,uguess,dt)
global K1 K2 Z  


delta = 1e-8;

for ii=1:Z+K1+2
   
    utemp = uguess;
    utemp(ii) = uguess(ii) + delta;
    resid =getresid(uold,uguess,dt);
    residDelta = getresid(uold,utemp,dt);
    jacob(:,ii) = (residDelta-resid)/delta;
    
end
%jacob(Z+K1+3,Z+K1+3) = 1;