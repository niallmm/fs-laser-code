function [ResIRFmovT] = getResidIRFmovT(oldh,hguess,dt,phi,TintOld,TintGuess,CintLold,CintLguess,CintSold,CintSguess, stateSold,stateSguess)% TO DO CODING OF A verson
global Mu1nm
% INTERFACE RESPONSE FUNCTIONS
% (a) Temperature as a fct. of front velocity - in the documentation
% \tilde{f}
% function of hdot, interface Temp, concentrations a interface and the
% state (amorphous / crystaline) of the solid



% % Here we need the real data - for now a 'fake' version
% fakeparam = 1e-3;
% fakeparam2 = 0.0;
% % test1
% 
% ResIRFmovT = fakeparam*(hguess - oldh)/dt +  ...
%     ( phi*    (1 - fakeparam2*CintLguess - TintGuess )+ ...
%      (1-phi)* (1 - fakeparam2*CintLold   - TintOld)...
%     );

% Equilibrium - interface at melting temperature
%ResIRFmovT = TintGuess - 1.0;


%
% NON- Equilibirum
% 

paramM = 0.0;% > 0 the absolute value of the slope: Tmelt + 1 - m*C
% have to normalize paramMu as paramMu = Mu/(L/tau/theta) where 
% theta = Tmelt-Tamb
%paramMu = 92; %this was for crystaline silicon... from hoglund
% theta = 1385; % Tamb = 300;  Tmelt = 1685; K
% L = 10; % nm
% paramMu = 0.267*theta/L;
paramMu = 36.9795; % for amorphous
%paramMu = Mu1nm;
% paramMu = 2*36;
ResIRFmovT =  (hguess - oldh) ...
       - dt * paramMu * ...
       ( phi    * (TintGuess - 1 + paramM * CintLguess) ...
       +(1-phi) * (TintOld   - 1 + paramM * CintLold  )    );   