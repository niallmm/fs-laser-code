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


kappaE = 1e-4;
vd  = vdglobal; % in 10nm/ns = 10 m/s

global check_IRFC

if (check_IRFC)
% option 1 
 % %CintLguess
 JIRFmovC(1, Z+K1) =   (hguess - oldh) * (CintSold - CintLold) ...
         -dt*vd*(   phi  * kappaE * (CintSold - CintLold) ...
               + (1-phi) * (CintSold   - kappaE * CintLold)    );                            
 % %CintSguess
 JIRFmovC(1, Z+K1+1)= -(hguess - oldh) * (CintSold - CintLold) ...
         +dt*vd*(   phi           * (CintSold - CintLold) ...
               + (1-phi) * (CintSold   - kappaE * CintLold)    );                            
 %      
 % % hguess
 JIRFmovC(1, Z+K1+2) = - (CintSguess - CintLguess) * (CintSold - CintLold);
 
else
% option 2    
 %CintLguess
 JIRFmovC(1, Z+K1)   = - dt * vd * kappaE + (hguess - oldh);

 %CintSguess
 JIRFmovC(1, Z+K1+1) =   dt * vd          - (hguess - oldh);

 % hguess
 JIRFmovC(1, Z+K1+2) = - (CintSguess - CintLguess)* kappaE;
 
end
