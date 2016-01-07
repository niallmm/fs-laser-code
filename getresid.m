function [resid] = getresid(uold,uguess,dt)
global  K1 K2 Z xmesh dx1 dx2 phi didi_1 dwa dwb Aparam

global TguessL TguessS CguessL CguessS ToldL ToldS ColdL ColdS tglobal 
global solidStateGlobal % concSolidEnforceGlobal

global TintFct_global TintInvLam_global K2as_NND time_NND


hguess = uguess(end);
oldh = uold(end);


TguessL = uguess(1:K1);
TguessS = uguess(K1+1:K1+K2);
CguessL = uguess(Z+1:Z+K1);
CguessS = uguess(Z+K1+1);

ToldL = uold(1:K1);
ToldS = uold(K1+1:K1+K2);
ColdL = uold(Z+1:Z+K1);
ColdS = uold(Z+K1+1);



% need to determine the concentration in the solid phase as well as he
% status (amorpheous / crysalline) from interpolated vector in the lab
% frame 



% CHECK: at the moment CguessSfull and stateSguessfull and thereby the new
% Diffusion constants depend on hguess. Since the derivative of all these
% parameters with respect to hguess is not analytically available, the
% Jacobian d(.)/dhguess is missing terms.
% Solution: We assume that the state and the concentration in the solid
% vary smoothly and use the OLD values to compute the diffusion constants:


%ColdSfull   = getCSfromCSlab(oldh,ColdS); % needs also the various meshes,.... - provided via global variables
% Value at x=h taken from the actual state vector
%ColdSfull(1) = ColdS;

% CHECK
%CguessSfull  = getCSfromCSlab(hguess,CguessS); % needs also the various meshes,.... - provided via global variables
%CguessSfull = ColdSfull;

% for hdot>0 (melting) this sets the conc. at the moving interface inside
% the solid (see ResICmovC)

% We enforce this solid concentration when hdot >0
% CSenforce = concSolidEnforceGlobal;

% Value at x=h taken from the actual state vector
%CguessSfull(1) = CguessS;

% The state

% NOTE: set both states to the globally enforced value... this avoids
% having to deal with non-differential variations of the diffusion
% constants with the state... 
stateSoldfull     = getStateSfromstateLab(oldh);
stateSoldfull(1)  = solidStateGlobal;

% CHECK
%stateSguessfull   = getStateSfromstateLab(hguess); % also requires the various meshes and the quantities in lab space - provided via global variables 
stateSguessfull   = stateSoldfull;
stateSguessfull(1)= solidStateGlobal; 


% at the interface
stateSold   = stateSoldfull(1);
stateSguess = stateSguessfull(1);


% Now compute the thermal diffusion constant everywhere in the solid

% NO DEPENDENCE ON CONCENTRATION!!!
%DToldS   = getDT(ToldS,ColdSfull,stateSoldfull);
%DTguessS = getDT(TguessS,CguessSfull,stateSguessfull);

DToldS   = getDTS(ToldS,stateSoldfull);
DTguessS = getDTS(TguessS,stateSguessfull);

% Thermal diffusion constant in the liquid
DToldL   = getDTL(ToldL,ColdL);
DTguessL = getDTL(TguessL,CguessL);

% Solute diffusion constants in the liquid phase

DColdL   = getDCL(ToldL,ColdL);
DCguessL = getDCL(TguessL,CguessL);



% Thermal diffusivity depending on interface temperature,.... at the
% interface - index K1, K1+1
KappaThermguessL = getKappaTherm(TguessL(K1),0); % State = 0 liquid
KappaThermguessS = getKappaTherm(TguessS(1),stateSguess);
KappaThermoldL   = getKappaTherm(ToldL(K1)  ,0); % State = 0 liquid
KappaThermoldS   = getKappaTherm(ToldS(1),  stateSold);




% Diffusion constant for the solute in the liquid phase - at the interface

DSguess = DCguessL(K1);
DSold   = DColdL(K1);

% flag that ckecks if we are melting (hdot >0)
hincrease = ((uold(K1) - uold(K1-1))*KappaThermoldL/dx1 - (uold(K1+2) - uold(K1+1))*KappaThermoldS/dx2(1)) <0;



if (hincrease)

% Residul for temperature in the liquid phase (uniform grid)


ResTliq = dx1*(TguessL(2:K1-1) - ToldL(2:K1-1)) - (hguess-oldh) * xmesh(2:K1-1) .*...
    ( phi/hguess   * (TguessL(3:K1)-TguessL(2:K1-1)) ...
    + (1-phi)/oldh * (ToldL(3:K1)  -ToldL(2:K1-1)  )...
    ) ...
    - dt./(2.0 * dx1) .* ...
    ( phi/(hguess*hguess)   * (  TguessL(3:K1)   .* ( DTguessL(3:K1) + DTguessL(2:K1-1) ) ...
                               - TguessL(2:K1-1) .* ( DTguessL(3:K1) + 2.0 * DTguessL(2:K1-1) + DTguessL(1:K1-2) )  ...
                               + TguessL(1:K1-2) .* ( DTguessL(2:K1-1) + DTguessL(1:K1-2) ) ...
                              ) ...
     +(1-phi) / (oldh*oldh) * (  ToldL(3:K1)     .* ( DToldL(3:K1)   + DToldL(2:K1-1)   ) ...
                               - ToldL(2:K1-1)   .* ( DToldL(3:K1)   + 2.0 * DToldL(2:K1-1)   + DToldL(1:K1-2)   )  ...
                               + ToldL(1:K1-2)   .* ( DToldL(2:K1-1)   + DToldL(1:K1-2)   ) ...
                              ) ...
    );                         
    


% Residul for temperature in the solid phase (non-uniform grid)
% NOTE: MAYBE DEVIDE BY ./xmesh(K1+2:K1+K2-1) 


% 

%K2as = getK2as_NND(oldh,tglobal);
%K2as_NND = K2as
K2as = K2as_NND;

%largeparam = 10;
%K2as = find(oldh*xmesh(K1+1:K1+K2)- largeparam * sqrt(4*Dconst_BBC*tglobal,1,'first');


% for 2:K2as-1: real residual
% for K2as:K2-1: set residual based on initial temperature

xeval = hguess .* xmesh(K1+K2as:K1+K2-1);
%xeval = oldh .* xmesh(K1+K2as:K1+K2-1);

one_over_lambda = TintInvLam_global(xeval);

InitTemp = TintFct_global(xeval);

% Ttarget = InitTemp ...
%                   * exp(DTguessS(2:K2as-1) .* one_over_lambda.^2 .*tglobal);
%  scaleparam = 1e-2;   
%  ResTsolAs = scaleparam * (TguessS(K2as:K2-1) - InitTemp ...
%                    .* exp(DTguessS(K2as:K2-1) .* one_over_lambda.^2 .*tglobal));
  time = time_NND;  
 % time = tglobal;  
  ResTsolAs = TguessS(K2as:K2-1) - InitTemp ...
                    .* exp(DTguessS(K2as:K2-1) .* one_over_lambda.^2 .*time);
   
               %global Dconst_BBC
%ResTsolAs = TguessS(K2as:K2-1) - InitTemp ...
 %                 .* exp(DToldS(K2as:K2-1) .* one_over_lambda.^2 .*tglobal);

ResTsolFull = dx2(2:K2as-1) .*(TguessS(2:K2as-1) - ToldS(2:K2as-1)) ...
     - (hguess-oldh) .* xmesh(K1+2:K1+K2as-1) .*...
     ( phi/hguess      * ( TguessS(3:K2as)- TguessS(2:K2as-1) )  ...
       + (1-phi)/oldh  * ( ToldS(3:K2as)  - ToldS(2:K2as-1)   )  ...
     ) ...
     - dt./( dx2(2:K2as-1) + dx2(1:K2as-2)) .* ...
     ( phi/(hguess*hguess)   * (  TguessS(3:K2as)   .* ( DTguessS(3:K2as) + DTguessS(2:K2as-1) ) ...
                                - TguessS(2:K2as-1) .* ( DTguessS(3:K2as) + DTguessS(2:K2as-1) + didi_1(2:K2as-1) .* ( DTguessS(2:K2as-1) + DTguessS(1:K2as-2) ) ) ...
                                + TguessS(1:K2as-2) .*                                       didi_1(2:K2as-1) .* ( DTguessS(2:K2as-1) + DTguessS(1:K2as-2) ) ...
                               ) ...
      +(1-phi) / (oldh*oldh) * (  ToldS(3:K2as)     .* ( DToldS(3:K2as)   + DToldS(2:K2as-1)   ) ...
                                - ToldS(2:K2as-1)   .* ( DToldS(3:K2as)   + DToldS(2:K2as-1)   + didi_1(2:K2as-1) .* ( DToldS(2:K2as-1)   + DToldS(1:K2as-2)   ) ) ...
                                + ToldS(1:K2as-2)   .*                                       didi_1(2:K2as-1) .* ( DToldS(2:K2as-1)   + DToldS(1:K2as-2)   ) ...
                               ) ...
     );                         

              

ResTsol = [ResTsolFull;ResTsolAs];

% residtest =    ResTsolFull(end)  
%    residtest2 =    ResTsolAs(1)  
   
% ResTsol = dx2(2:K2-1) .*(TguessS(2:K2-1) - ToldS(2:K2-1)) ...
%     - (hguess-oldh) .* xmesh(K1+2:K1+K2-1) .*...
%     ( phi/hguess      * ( TguessS(3:K2)- TguessS(2:K2-1) )  ...
%       + (1-phi)/oldh  * ( ToldS(3:K2)  - ToldS(2:K2-1)   )  ...
%     ) ...
%     - dt./( dx2(2:K2-1) + dx2(1:K2-2)) .* ...
%     ( phi/(hguess*hguess)   * (  TguessS(3:K2)   .* ( DTguessS(3:K2) + DTguessS(2:K2-1) ) ...
%                                - TguessS(2:K2-1) .* ( DTguessS(3:K2) + DTguessS(2:K2-1) + didi_1(2:K2-1) .* ( DTguessS(2:K2-1) + DTguessS(1:K2-2) ) ) ...
%                                + TguessS(1:K2-2) .*                                       didi_1(2:K2-1) .* ( DTguessS(2:K2-1) + DTguessS(1:K2-2) ) ...
%                               ) ...
%      +(1-phi) / (oldh*oldh) * (  ToldS(3:K2)     .* ( DToldS(3:K2)   + DToldS(2:K2-1)   ) ...
%                                - ToldS(2:K2-1)   .* ( DToldS(3:K2)   + DToldS(2:K2-1)   + didi_1(2:K2-1) .* ( DToldS(2:K2-1)   + DToldS(1:K2-2)   ) ) ...
%                                + ToldS(1:K2-2)   .*                                       didi_1(2:K2-1) .* ( DToldS(2:K2-1)   + DToldS(1:K2-2)   ) ...
%                               ) ...
%     );                         





% Residul for the solute concentration in the liquid phase (uniform grid)


ResCliq = dx1*(CguessL(2:K1-1) - ColdL(2:K1-1)) - (hguess-oldh) * ...
    ( phi/hguess   * xmesh(2:K1-1).*(CguessL(3:K1)-CguessL(2:K1-1)) ...
    + (1-phi)/oldh * xmesh(2:K1-1) .* (ColdL(3:K1)  -ColdL(2:K1-1))...
    ) ...
    - dt./(2.0 * dx1) .* ...
    ( phi/(hguess*hguess)   * (  CguessL(3:K1)   .* ( DCguessL(3:K1) + DCguessL(2:K1-1) ) ...
                               - CguessL(2:K1-1) .* ( DCguessL(3:K1) + 2.0 * DCguessL(2:K1-1) + DCguessL(1:K1-2) )  ...
                               + CguessL(1:K1-2) .* ( DCguessL(2:K1-1) + DCguessL(1:K1-2) ) ...
                              ) ...
     +(1-phi) / (oldh*oldh) * (  ColdL(3:K1)     .* ( DColdL(3:K1)   + DColdL(2:K1-1)   ) ...
                               - ColdL(2:K1-1)   .* ( DColdL(3:K1)   + 2.0 * DColdL(2:K1-1)   + DColdL(1:K1-2)   )  ...
                               + ColdL(1:K1-2)   .* ( DColdL(2:K1-1)   + DColdL(1:K1-2)   ) ...
                              ) ...
    );                         




else 
% Resolidification - hdot <0

% Residul for temperature in the liquid phase (uniform grid)


ResTliq = dx1*(TguessL(2:K1-1) - ToldL(2:K1-1)) - (hguess-oldh) * xmesh(2:K1-1).*...
    ( phi/hguess   * (TguessL(2:K1-1)-TguessL(1:K1-2)) ...
    + (1-phi)/oldh * (ToldL(2:K1-1)  -ToldL(1:K1-2))...
    ) ...
    - dt./(2.0 * dx1) .* ...
    ( phi/(hguess*hguess)   * (  TguessL(3:K1)   .* ( DTguessL(3:K1) + DTguessL(2:K1-1) ) ...
                               - TguessL(2:K1-1) .* ( DTguessL(3:K1) + 2.0 * DTguessL(2:K1-1) + DTguessL(1:K1-2) )  ...
                               + TguessL(1:K1-2) .* ( DTguessL(2:K1-1) + DTguessL(1:K1-2) ) ...
                              ) ...
     +(1-phi) / (oldh*oldh) * (  ToldL(3:K1)     .* ( DToldL(3:K1)   + DToldL(2:K1-1)   ) ...
                               - ToldL(2:K1-1)   .* ( DToldL(3:K1)   + 2.0 * DToldL(2:K1-1)   + DToldL(1:K1-2)   )  ...
                               + ToldL(1:K1-2)   .* ( DToldL(2:K1-1)   + DToldL(1:K1-2)   ) ...
                              ) ...
    );                         




% Residul for temperature in the solid phase (non-uniform grid)
% NOTE: MAYBE DEVIDE BY ./xmesh(K1+2:K1+K2-1) 


ResTsol = dx2(2:K2-1) .*(TguessS(2:K2-1) - ToldS(2:K2-1)) ...
    - (hguess-oldh) * didi_1(2:K2-1) .* xmesh(K1+2:K1+K2-1) .*...
    ( phi/hguess      * ( TguessS(2:K2-1)- TguessS(1:K2-2) )  ...
      + (1-phi)/oldh  * ( ToldS(2:K2-1)  - ToldS(1:K2-2)   )  ...
    ) ...
    - dt./( dx2(2:K2-1) + dx2(1:K2-2)) .* ...
    ( phi/(hguess*hguess)   * (  TguessS(3:K2)   .* ( DTguessS(3:K2) + DTguessS(2:K2-1) ) ...
                               - TguessS(2:K2-1) .* ( DTguessS(3:K2) + DTguessS(2:K2-1) + didi_1(2:K2-1) .* ( DTguessS(2:K2-1) + DTguessS(1:K2-2) ) ) ...
                               + TguessS(1:K2-2) .*                                       didi_1(2:K2-1) .* ( DTguessS(2:K2-1) + DTguessS(1:K2-2) ) ...
                              ) ...
     +(1-phi) / (oldh*oldh) * (  ToldS(3:K2)     .* ( DToldS(3:K2)   + DToldS(2:K2-1)   ) ...
                               - ToldS(2:K2-1)   .* ( DToldS(3:K2)   + DToldS(2:K2-1)   + didi_1(2:K2-1) .* ( DToldS(2:K2-1)   + DToldS(1:K2-2)   ) ) ...
                               + ToldS(1:K2-2)   .*                                       didi_1(2:K2-1) .* ( DToldS(2:K2-1)   + DToldS(1:K2-2)   ) ...
                              ) ...
    );                         





% Residul for the solute concentration in the liquid phase (uniform grid)


ResCliq = dx1*(CguessL(2:K1-1) - ColdL(2:K1-1)) - (hguess-oldh) * xmesh(2:K1-1).*...
    ( phi/hguess   * (CguessL(2:K1-1)-CguessL(1:K1-2)) ...
    + (1-phi)/oldh * (ColdL(2:K1-1)  -ColdL(1:K1-2)  )...
    ) ...
    - dt./(2.0 * dx1) .* ...
    ( phi/(hguess*hguess)   * (  CguessL(3:K1)   .* ( DCguessL(3:K1) + DCguessL(2:K1-1) ) ...
                               - CguessL(2:K1-1) .* ( DCguessL(3:K1) + 2.0 * DCguessL(2:K1-1) + DCguessL(1:K1-2) )  ...
                               + CguessL(1:K1-2) .* ( DCguessL(2:K1-1) + DCguessL(1:K1-2) ) ...
                              ) ...
     +(1-phi) / (oldh*oldh) * (  ColdL(3:K1)     .* ( DColdL(3:K1)   + DColdL(2:K1-1)   ) ...
                               - ColdL(2:K1-1)   .* ( DColdL(3:K1)   + 2.0 * DColdL(2:K1-1)   + DColdL(1:K1-2)   )  ...
                               + ColdL(1:K1-2)   .* ( DColdL(2:K1-1)   + DColdL(1:K1-2)   ) ...
                              ) ...
    );                         


    
end



% Interface conditions at the moving front
% (a) energy conservation



ResBCmovA = hguess - oldh + ...
        dt*Aparam*phi/hguess   * (  KappaThermguessL * (1.5 * TguessL(K1) - 2.0     * TguessL(K1-1) + 0.5 * TguessL(K1-2) )/dx1 ...
                                  - KappaThermguessS * (-dwb * TguessS(3) +(dwa+dwb)* TguessS(2)    - dwa *   TguessS(1)  )      ) + ...
        dt*Aparam*(1-phi)/oldh * (  KappaThermoldL *   (1.5 * ToldL(K1) - 2.0       * ToldL(K1-1)   + 0.5 * ToldL(K1-2)   )/dx1 ...
                                  - KappaThermoldS *   (-dwb * ToldS(3) +(dwa+dwb)* ToldS(2)    - dwa *   ToldS(1)  )      );



% (b) solute conservation 



ResBCmovB = 0.5 * dx1 * ( hguess*hguess - oldh*oldh ) * ...
            (phi* (CguessL(K1) - CguessS(1) )  + (1-phi)*(ColdL(K1) - ColdS(1) )) ...
          + dt * (  DSguess * phi  * ( 1.5 * CguessL(K1) - 2.0 * CguessL(K1-1) + 0.5 * CguessL(K1-2) ) ...
                 +  DSold * (1-phi)* ( 1.5 * ColdL(K1)   - 2.0 * ColdL(K1-1)   + 0.5 * ColdL(K1-2)   )  );


% (c) continuity of T - index K1+1

ResBCmovC = TguessL(K1) - TguessS(1);


      

% INTERFACE RESPONSE FUNCTIONS
% (a) Temperature as a fct. of front velocity - in the documentation
% \tilde{f}
% function of hdot, interface Temp, concentrations a interface and the
% state (amorphous / crystaline) of the solid
ResIRFmovT = getResidIRFmovT(oldh,hguess,dt,phi,ToldL(K1),TguessL(K1),ColdL(K1),CguessL(K1),ColdS(1),CguessS(1), stateSold,stateSguess);% TO DO CODING OF A verson


% ODE for h in case of a non-equilibrium model
% for the pure equilibrium case it turns into a simple condition on thw
% interface temperature

% (b) jump in the concentration - in the documentation
% \tilde{g}  (only relevant for resolidification part - \dot{h} < 0
% function of hdot, interface Temp, concentrations a interface and the
% state (amorphous / crystaline) of the solid
if(~hincrease)
ResIRFmovC = getResidIRFmovC(oldh,hguess,dt,phi,ToldL(K1),TguessL(K1),ColdL(K1),CguessL(K1),ColdS(1),CguessS(1), stateSold,stateSguess);% TO DO CODING OF A verson


% (b2) for \dot{h} >0 (melting) we instead keep the solid concentration
% value at the moving interface at whatever value of clab(h) -
% interpolation


else
CSenforce = getInterfaceCSfromCSlab(hguess); 
ResICmovC = CguessS(1) - CSenforce;
end


% Boundary conditions at surface (x=0)
% (a) Temperature at the surface - in the documentation \tilde{F} (capital
% F)
% depends on TguessL(2), TguessL(1), hguess, tglobal (explicit time
% dependence 
ResBCsurfT = getResidBCsurfT(TguessL(1),TguessL(2),hguess,tglobal);



% (b) Concentration at the surface - in the documentation \tilde{G} (capital G)
% NOTE: Can be explicitly time dependent - provide time as a parameter -
% not as an element of the state vector

% depends on CguessL(2), CguessL(1), Tguess(1) - surface temperature,
% hguess, tglobal
ResBCsurfC = getResidBCsurfC(CguessL(1),CguessL(2),TguessL(1),hguess,tglobal);

         

% Boundary conditions in the bulk x \to \infty (actually L2)
% Temperature approaches bulk temperature (=0 in our units). Alternatively
% we can consider an effective BC to account fo the finite value of L2 - in
% the documentation \tilde{H}
% NOTE: MAYBE DEVELOP THAT EFFECTIVE BC

%ResBCbulkT = TguessS(K2);


if (hincrease)
ResBCbulkT = getResidBCbulkT(TguessS(K2),hguess,tglobal);
else
ResBCbulkT = getResidBCbulkTresold(dt,oldh,hguess,TguessS,tglobal,ToldS(K2-3));    
end

% Put everything together
% Construct the residual - the order is somewhat historically determined
% and loosely reflects which element of the statevector is fixed by setting
% the resiudal to zero


if (hincrease)
resid = [ResBCsurfT; ...    % temp at surface
        ResTliq; ...        % temp diffusion in the liquid
        ResBCmovA; ...      % energy conservation at the moving front
        ResBCmovC; ...      % continuity of the temperature across the interface
        ResTsol; ...        % temp diffusion in the solid
        ResBCbulkT; ...     % temp at x \to \ifty
        ResBCsurfC; ...     % solute concentration at the surface 
        ResCliq; ...        % solute diffusion in the liquid
        ResBCmovB; ...      % solute conservation at the moving front
        ResICmovC; ...      % solute conc. at the front (in the solid). From initial conc. profile - only for \dot{h}>0
        ResIRFmovT];%...      % Interface reponse function - basically interface temperature vs. front velocity
else

resid = [ResBCsurfT; ...    % temp at surface
        ResTliq; ...        % temp diffusion in the liquid
        ResBCmovA; ...      % energy conservation at the moving front
        ResBCmovC; ...      % continuity of the temperature across the interface
        ResTsol; ...        % temp diffusion in the solid
        ResBCbulkT; ...     % temp at x \to \ifty
        ResBCsurfC; ...     % solute concentration at the surface 
        ResCliq; ...        % solute diffusion in the liquid
        ResBCmovB; ...      % solute conservation at the moving front
%        ResICmovC; ...      % solute conc. at the front (in the solid). From initial conc. profile - only for \dot{h}>0
        ResIRFmovC; ...     % Interface response function for the concentration jump at the moving front - only for \dot{h}<0
        ResIRFmovT];% ...        % Interface response function - basically interface temperature vs. front velocity

% JUST FOR CHECKING

%ResBCbulkT_CHECK = getResidBCbulkTresold(dt,oldh,hguess,TguessS,tglobal,ToldS(K2-3));
%error('stop in getresid')


end    
    

    

end