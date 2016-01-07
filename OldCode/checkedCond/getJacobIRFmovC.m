function [JIRFmovC] = getJacobIRFmovC(oldh,hguess,dt,phi,TintOld,TintGuess,CintLold,CintLguess,CintSold,CintSguess, stateSold,stateSguess)% TO DO CODING OF A verson
global K1  K2 Z
global kappaEglobal vdglobal % only for varying those parameters

% JACOBIAN - 
% Vector: T[1:K1] liquid, T[K1+1,K1+K2] solid, C[1:K1] liquid, C[K1+1] (SINGLE
% POINT) solid, h

% INTERFACE RESPONSE FUNCTIONS
% (b) jump in the concentration - in the documentation
% \tilde{g}  (only relevant for resolidification part - \dot{h} < 0
% function of hdot, interface Temp, concentrations a interface and the
% state (amorphous / crystaline) of the solid

JIRFmovC = sparse(1,Z+K1+2);

      
%fakeparam = 1e0;

% % Interface temperature - index K1
% %JIRFmovC(1,K1) = - phi*dt ;
% % Interface concen in the liquid - index Z+K1
% JIRFmovC(1,Z+K1) = (oldh-hguess) ;
% % Interface concen in the solid - index Z+K1+1
% JIRFmovC(1,Z+K1+1) = -fakeparam;%- sqrt(fakeparam) ;
% % h guess - index Z+K1+2
% JIRFmovC(1,Z+K1+2) = -CintLguess;
% 
% Interface conc in the liquid
%JIRFmovC(1,Z+K1) = 1.0;
% Interface conc in the sold
%JIRFmovC(1,Z+K1+1) = -fakeparam;
% hguess
%JIRFmovC(1,Z+K1+2) = fakeparam*CintSguess;

% =======================================
% Equilibrium Cond.
% =======================================

% kappa = 1; %find papers k(hdot, CL, CS) 
% 
% % Interface concentration (last point) in the liquid - index Z+K1
% JIRFmovC(1, Z+K1) = 1/sqrt(kappa);
% % Interface concentration (first point) in the solid - index Z+K1+1
% % (TMS) should be a minus sign
% JIRFmovC(1, Z+K1+1) =  - sqrt(kappa);
% 

% =======================================
% Non Equilibrium Cond.
% =======================================

% kappaE = 1e-4;
% %CintLguess
% JIRFmovC(1, Z+K1) =   (hguess - oldh) * (CintSold - CintLold) ...
%         - dt * (   phi  * kappaE * (CintSold - CintLold) ...
%               + (1-phi) * (CintSold   - kappaE * CintLold)    );                            
% %CintSguess
% JIRFmovC(1, Z+K1+1)= -(hguess - oldh) * (CintSold - CintLold) ...
%         + dt * (   phi           * (CintSold - CintLold) ...
%               + (1-phi) * (CintSold   - kappaE * CintLold)    );                            
%      
% % hguess
% JIRFmovC(1, Z+K1+2) = - (CintSguess - CintLguess) * (CintSold - CintLold);
                               
% =======================================
% Non - Equilibrium Cond.
% new discretization (06/26/11) because the old one generates a row of the
% jacobian which identically vanishes for vanishing concentrations. This
% causes the jacobian to be singular!!!!
% Here: simpler discritization - potentially unstable wrt to time
% integration - let's hope ;-)
% =======================================
% note: in our units V_D = 1
                                
%  kappaE = 1e-4;
%  vd     = 1.0; % in 10nm/ns = 10 m/s
% Hacking - define them via global variables
kappaE = kappaEglobal;
vd     = vdglobal; % in 10nm/ns = 10 m/s


%CintLguess
JIRFmovC(1, Z+K1)   = - dt * vd * kappaE + (hguess - oldh);

%CintSguess
JIRFmovC(1, Z+K1+1) =   dt * vd          - (hguess - oldh);

% hguess
JIRFmovC(1, Z+K1+2) = - (CintSguess - CintLguess);
