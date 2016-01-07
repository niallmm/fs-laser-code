function [ResIRFmovC] = getResidIRFmovC(oldh,hguess,dt,phi,TintOld,TintGuess,CintLold,CintLguess,CintSold,CintSguess, stateSold,stateSguess)% TO DO CODING OF A verson
global kappaEglobal vdglobal % only for varying those parameters

% INTERFACE RESPONSE FUNCTIONS

% (b) jump in the concentration - in the documentation
% \tilde{g}  (only relevant for resolidification part - \dot{h} < 0
% function of hdot, interface Temp, concentrations a interface and the
% state (amorphous / crystaline) of the solid




Cthreshold = 1e-6;
kappaE = 1e-4;
vd = vdglobal; % in 10nm/ns = 10 m/s

global check_IRFC
check_IRFC = CintLold > Cthreshold;
% note  make this a global variables so that the jacobian is evaluated for
% exactly the same option


if (check_IRFC)
% option 1
  ResIRFmovC = - (hguess - oldh) * (CintSguess - CintLguess) ...
                              * (CintSold - CintLold) ...
     + dt*vd*(   phi    * (CintSguess - kappaE * CintLguess) ...
                                     * (CintSold - CintLold) ...
              + (1-phi) * (CintSold   - kappaE * CintLold)  ...
                                     * (CintSguess - CintLguess) ...
             );  
else
%option 2    
  ResIRFmovC =  dt * vd * (CintSguess - kappaE * CintLguess) ...
                        - (hguess - oldh) * (CintSguess - CintLguess) ;  
    
    
end


 return            