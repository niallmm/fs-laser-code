function [ResIRFmovC] = getResidIRFmovC(oldh,hguess,dt,phi,TintOld,TintGuess,CintLold,CintLguess,CintSold,CintSguess, stateSold,stateSguess)% TO DO CODING OF A verson
global kappaEglobal vdglobal % only for varying those parameters

% INTERFACE RESPONSE FUNCTIONS

% (b) jump in the concentration - in the documentation
% \tilde{g}  (only relevant for resolidification part - \dot{h} < 0
% function of hdot, interface Temp, concentrations a interface and the
% state (amorphous / crystaline) of the solid


% Here we need the real data - for now a 'fake' version
%fakeparam = 1e0;
%ResIRFmovC = CintLguess./sqrt(fakeparam) - sqrt(fakeparam)*CintSguess;
%
%ResIRFmovC = (oldh-hguess)*CintLguess - fakeparam*CintSguess;

%ResIRFmovC = 0;%CintSguess;

% =======================================
% Equilibrium Cond.
% =======================================

%  kappa = 1; %find papers k(hdot, CL, CS) 
% % 
%  ResIRFmovC = CintLguess./sqrt(kappa) - sqrt(kappa)*CintSguess;


% =======================================
% Non - Equilibrium Cond.
% =======================================
% note: in our units V_D = 1

% kappaE = 1e-4;
% ResIRFmovC = - (hguess - oldh) * (CintSguess - CintLguess) ...
%                              * (CintSold - CintLold) ...
%     + dt * (   phi    * (CintSguess - kappaE * CintLguess) ...
%                                     * (CintSold - CintLold) ...
%             + (1-phi) * (CintSold   - kappaE * CintLold)  ...
%                                     * (CintSguess - CintLguess) ...
%            );
%

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
 
 
 
 ResIRFmovC = dt * vd * (CintSguess - kappaE * CintLguess) ...
                      - (hguess - oldh) * (CintSguess - CintLguess) ;
             
             
             