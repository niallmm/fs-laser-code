function [JBCsurfC] = getJacobBCsurfC(CguessL1, CguessL2, Tsurface,hguess, tglobal)
global K1  K2 Z cflux dx1
% Vector: T[1:K1] liquid, T[K1+1,K1+K2] solid, C[1:K1] liquid, C[K1+1] (SINGLE
% POINT) solid, h


% (b) Concentration at the surface - in the documentation \tilde{G} (capital G)
% NOTE: Can be explicitly time dependent - provide time as a parameter -
% not as an element of the state vector
% depends on CguessL(2), CguessL(1), Tguess(1) - surface temperature,
% hguess, tglobal
JBCsurfC = sparse(1,Z+K1+2);
% for const at x-0
JBCsurfC(1,Z+1) = -1.;
JBCsurfC(1,Z+2) = 1.;
% for constant flux at x=0
% 
% JBCsurfC = sparse(1,Z+K1+2);
sulfFlux= cflux;
  Dliquid  = getDCL(1,0);
JBCsurfC(1,Z+1) = -Dliquid;
JBCsurfC(1,Z+2) = Dliquid;
JBCsurfC(1,Z+K1+2) = sulfFlux*dx1;
% 
% % CHECK
% % Tsurface
% 
 % JBCsurfC(1,1) =  CguessL1;
% % CguessL1
 % JBCsurfC(1,Z+1) = Tsurface + hguess;
% % CguessL2
% % 
 % JBCsurfC(1,Z+2) = - hguess;
% % 
% % hguess
% JBCsurfC(1,Z+K1+2) = -(CguessL2-CguessL1);



