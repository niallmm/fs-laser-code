function [ResIRFmovT] = getResidIRFmovT(oldh,hguess,dt,phi,TintOld,TintGuess,CintLold,CintLguess,CintSold,CintSguess, stateSold,stateSguess)% TO DO CODING OF A verson
% global Mu1nm
% INTERFACE RESPONSE FUNCTIONS
% (a) Temperature as a fct. of front velocity - in the documentation
% \tilde{f}
% function of hdot, interface Temp, concentrations a interface and the
% state (amorphous / crystaline) of the solid


%
% NON- Equilibirum
% 

paramM = 0.0;% > 0 the absolute value of the slope: Tmelt + 1 - m*C
% this parameter would make the inteface velocity dependent on the solute
% conentration. The simplist version is that it doesn't

% from Hogland mu = 0.267 m/(s K) for amorphous 
% and mu = 0.0667 m/(s K) for crystaline
% have to normalize paramMu as paramMu = Mu/(L/tau/theta) where 
% theta = Tmelt-Tamb = 1385;
% L = 10; % nm
paramMu = 92; %this was for crystaline silicon... from hoglund
%paramMu = 37.0; % for amorphous

ResIRFmovT =  (hguess - oldh) ...
       - dt * paramMu * ...
       ( phi    * (TintGuess - 1 + paramM * CintLguess) ...
       +(1-phi) * (TintOld   - 1 + paramM * CintLold  )    );   