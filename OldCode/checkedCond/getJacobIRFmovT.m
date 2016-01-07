function [JIRFmovT] = getJacobIRFmovT(oldh,hguess,dt,phi,TintOld,TintGuess,CintLold,CintLguess,CintSold,CintSguess, stateSold,stateSguess)% TO DO CODING OF A verson
global K1  K2 Z Mu1nm

% JACOBIAN - 
% Vector: T[1:K1] liquid, T[K1+1,K1+K2] solid, C[1:K1] liquid, C[K1+1] (SINGLE
% POINT) solid, h

% INTERFACE RESPONSE FUNCTIONS
% (a) Temperature as a fct. of front velocity - in the documentation
% \tilde{f}
% function of hdot, interface Temp, concentrations a interface and the
% state (amorphous / crystaline) of the solid

JIRFmovT = sparse(1,Z+K1+2);  

% %Temperature on l/sol boundary
% %JBCmovT(1,K1) = -1;
% %Concentration on l/sol boundary (in liquid)
% %JBCmovT(1,Z+K1) = - mphase;
% 
% % test1
% 
% fakeparam = 1e-3;
% fakeparam2 = 0.0;
% 
% % Interface temperature - index K1
% JIRFmovT(1,K1) = - phi ;
% % Interface concen in the liquid - index Z+K1
% JIRFmovT(1,Z+K1) = - phi*fakeparam2 ;
% % h guess - index Z+K1+2
% JIRFmovT(1,Z+K1+2) = fakeparam/dt ;


% equlilibirum
%JIRFmovT(1,K1) = 1.0 ;
%
% NON- Equilibirum
% 

paramM = 0.0;% > 0 the absolute value of the slope: Tmelt + 1 - m*C
% have to normalize paramMu as paramMu = Mu/(L/tau/theta) where 
% theta = Tmelt-Tamb
%paramMu = 9.2; %this was for crystaline silicon... from hoglund
% theta = 1385; % Tamb = 300;  Tmelt = 1685; K
% L = 10; % nm
% paramMu = 0.267*theta/L;
paramMu = 36.9795; % for amorphous .. from hoglund
% paramMu = 2*36;
%paramMu = Mu1nm;
% Interface temperature - index K1
JIRFmovT(1,K1)   = - dt *  paramMu * phi ;
% Interface concen in the liquid - index Z+K1
JIRFmovT(1,Z+K1) = - dt *  paramMu * phi * paramM;
% h guess - index Z+K1+2
JIRFmovT(1,Z+K1+2) =  1.0 ;