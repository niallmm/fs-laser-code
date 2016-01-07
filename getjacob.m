function [jacob] = getjacob(uold,uguess,dt)
global  K1 K2 Z xmesh dx1 dx2 phi didi_1 dwa dwb Aparam

global TguessL TguessS CguessL CguessS ToldL ToldS ColdL ColdS tglobal 

global solidStateGlobal 


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


%ColdSfull   = getCSfromCSlab(oldh,ColdS); % needs also the various meshes,.... - provided via global variables
% Value at x=h taken from the actual state vector
%ColdSfull(1) = ColdS;

% CHECK
%CguessSfull  = getCSfromCSlab(hguess,CguessS); % needs also the various meshes,.... - provided via global variables
%CguessSfull = ColdSfull;

%CguessSfull(1) = CguessS;


% The state

stateSoldfull     = getStateSfromstateLab(oldh);
stateSoldfull(1)  = solidStateGlobal;
% CHECK
%stateSguessfull   = getStateSfromstateLab(hguess); % also requires the various meshes and the quantities in lab space - provided via global variables 
stateSguessfull = stateSoldfull;

stateSguessfull(1)= solidStateGlobal; 

% at the interface
stateSold   = stateSoldfull(1);
stateSguess = stateSguessfull(1);


% Now compute the thermal diffusion constant everywhere in the solid

%DToldS   = getDT(ToldS,ColdSfull,stateSoldfull);
% NO DEPENDENCE ON CONCENTRATION!!!
%DTguessS = getDT(TguessS,CguessSfull,stateSguessfull);
DTguessS = getDTS(TguessS,stateSguessfull);


% Thermal diffusion constant in the liquid
%DToldL   = getDTL(ToldL,ColdL);
DTguessL = getDTL(TguessL,CguessL);


% Solute diffusion constants in the liquid phase

%DColdL   = getDCL(ToldL,ColdL);
DCguessL = getDCL(TguessL,CguessL);



% Thermal diffusivity depending on interface temperature,.... at the
% interface - index K1, K1+1
KappaThermguessL = getKappaTherm(TguessL(K1),0); % State = 0 liquid
KappaThermguessS = getKappaTherm(TguessS(1),stateSguess);
KappaThermoldL   = getKappaTherm(ToldL(K1)  ,0); % State = 0 liquid
KappaThermoldS   = getKappaTherm(ToldS(1),  stateSold);




% Diffusion constant for the solute in the liquid phase - at the interface

DSguess = DCguessL(K1);
%DSold   = DColdL(K1);

% flag that ckecks if we are melting (hdot >0)
hincrease = ((uold(K1) - uold(K1-1))*KappaThermoldL/dx1 - (uold(K1+2) - uold(K1+1))*KappaThermoldS/dx2(1)) <0;

% Partial derivatives of the diffusion constantd

dDTdTguessL = getDTdTL(TguessL,CguessL);

dDTdCguessL = getDTdCL(TguessL,CguessL);

% NO DEPENDENCE ON CONCENTRATION!!!
%dDTdTguessS = getDTdT(TguessS,CguessSfull,stateSguessfull);
dDTdTguessS = getDTSdT(TguessS,stateSguessfull);
% THIS VANISHES NOW
%dDTdCguessS = getDTdC(TguessS(1),CguessS,stateSguess);

dDCdCguessL = getDCdCL(TguessL,CguessL);

dDCdTguessL = getDCdTL(TguessL,CguessL);

dKappaThermdTguessL = getKappaThermdT(TguessL(K1),0); % State = 0 liquid

dKappaThermdTguessS = getKappaThermdT(TguessS(1),stateSguess);

dKappaThermdCguessL = getKappaThermdC(TguessL(K1),0); % State = 0 liquid

dKappaThermdCguessS = getKappaThermdC(TguessS(1),stateSguess);


if (hincrease) 

    % ** ** **
    % JTT 
    
    % Temperatur (liquid) inner points: JTTliqIN (K1-2) x K1
 
    % (JTTliquid) Jacobian temperature residuals with respect to temperature - JTT in the liquid phase (uniform grid) 

    
    % diagonal
    diag =  dx1 ...
            + phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( DTguessL(3:K1) +  2.0 * DTguessL(2:K1-1) + DTguessL(1:K1-2)   ) ...
            + (1 - oldh/hguess) * phi * xmesh(2:K1-1);


    
    % superdiagonal
    supdia= - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( DTguessL(3:K1) +  DTguessL(2:K1-1) ) ...
            - (1 - oldh/hguess) * phi * xmesh(2:K1-1);
 
    
    
    
    % subdiagonal
    subdia = - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( DTguessL(2:K1-1) +  DTguessL(1:K1-2) ); 
    
    

    % ADD Coupling via the Diffusion constant
    
    
    % Derivative with respect to DTguess
    diag2  = - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( TguessL(3:K1) -  2.0 * TguessL(2:K1-1) + TguessL(1:K1-2)   );
    
    supdia2= - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( TguessL(3:K1) -  TguessL(2:K1-1) );
                
    subdia2= - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( -  TguessL(2:K1-1) +  TguessL(1:K1-2) ); 
             
    % multiply with partical derivatives of DT with respect to temperature
    % NOTE: Dj with j the COLUMN index -> for the diag / sub / sup diag we
    % have to multiply with
    % diag2 = diag2 .* (2:K1-1) (i,j) 
    % supdia = supdia .* (3:K1) with j = i+1 (index of the diag
    % runsover i !!!!
    % subdia = subdia .* (1:K1-2) with j = i-1 (index of the diag
    % runsover i !!!!
    
    diag   = diag   + diag2   .* dDTdTguessL(2:K1-1);
    supdia = supdia + supdia2 .* dDTdTguessL(3:K1);
    subdia = subdia + subdia2 .* dDTdTguessL(1:K1-2);
    

    %JTTliqIN = spdiags([subdia diag supdia], 0:2, K1-2, K1) ;
    JTTliqIN = myspdiag(subdia,diag,supdia);
    % first and last row are missing s. t. diagonal elements are at (i,i+1)

                
    
    % JTC vanishes no more due to coupling via the state dependent diffusion constant
    % dimension (K1-2) x K1 
    % Derivative with respect to DTguess multiplied with partical
    % derivatives of DT with respect to concentration
    
    diag   = diag2   .* dDTdCguessL(2:K1-1);
    supdia = supdia2 .* dDTdCguessL(3:K1);
    subdia = subdia2 .* dDTdCguessL(1:K1-2);
    
    
    JTCliqIN = myspdiag(subdia,diag,supdia);
    
    
% N E W
K2as = K2as_NND;

% 1. for the real residuals from 2:K2as-1 which depend on 1:K2as
% -> ITTsolINfull (K2as-2) x K2as

% 2. does not depend on the remaining ones - add sparse(K2as-2,K2-K2as)

% 3. remaining residuals K2as:K2-1 - they don't depend on 1:K2as-1
%   sparse(K2-1-K2as,K2as-1)

% 4. diagonal for K2as:K2-1:   speye(K2-1-K2as,K2-1-K2as)
    
    % ** ** **
    % JTT 
     
    % Temperatur (solid) inner points: JTTsolIN (K2-2) x K2
 
    % (JTTsolid) Jacobian temperature residuals with respect to temperature - JTT in the solid phase (non-uniform grid) 
    % maybe devide by ./xmesh(K1+2:K1+K2-1) (last row xmesh ->
    % ones(K2-2,1);
    
    % NOTE: IN THE OLD IMPLEMENTATION OF THE DIAGONAL THE DIFFUSIVE TERM WAS PROB WRONG -
    % not devided by xmesh !!!!!
      
    diag = [dx2(2:K2as-1) ...
            + phi*dt/(hguess*hguess) ./(dx2(2:K2as-1) + dx2(1:K2as-2)) .* ...
                    ( DTguessS(3:K2as) +  (1+didi_1(2:K2as-1)) .* DTguessS(2:K2as-1) + didi_1(2:K2as-1) .* DTguessS(1:K2as-2)   ) ...
            + (1 - oldh/hguess) * phi * xmesh(K1+2:K1+K2as-1) ...
           ;  ones(K2-K2as,1) ...
           ];
            
%     checkJacob3 = max(abs(dx2(2:K2as-1) ...
%             + phi*dt/(hguess*hguess) ./(dx2(2:K2as-1) + dx2(1:K2as-2)) .* ...
%                     ( DTguessS(3:K2as) +  (1+didi_1(2:K2as-1)) .* DTguessS(2:K2as-1) + didi_1(2:K2as-1) .* DTguessS(1:K2as-2)   ) ...
%             + (1 - oldh/hguess) * phi * xmesh(K1+2:K1+K2as-1)))
%     

    % superdiagonal
    supdia= [- phi*dt/(hguess*hguess) ./(dx2(2:K2as-1) + dx2(1:K2as-2)) .* ...
                    ( DTguessS(3:K2as) +  DTguessS(2:K2as-1) ) ...
            - (1 - oldh/hguess) * phi * xmesh(K1+2:K1+K2as-1)...
            ; sparse(K2-K2as,1) ...
            ];
 
    
    
    %subdiagonal
    subdia = [- phi*dt/(hguess*hguess) ./(dx2(2:K2as-1) + dx2(1:K2as-2)) .* ...
                    didi_1(2:K2as-1) .* ( DTguessS(2:K2as-1) +  DTguessS(1:K2as-2) )...
             ; sparse(K2-K2as,1) ...
             ];
    

       % ADD Coupling via the Diffusion constant
    
    xeval = hguess .* xmesh(K1+K2as:K1+K2-1);
    InitTemp = TintFct_global(xeval);
    one_over_lambda = TintInvLam_global(xeval);
    time = time_NND;
  %  time = tglobal;
    Ttarget =  InitTemp ...
                   .* exp(DTguessS(K2as:K2-1) .* one_over_lambda.^2 .*time); 

    % Derivative with respect to DTguess
    diag2  = [- phi*dt/(hguess*hguess) ./(dx2(2:K2as-1) + dx2(1:K2as-2)) .* ...
                    ( TguessS(3:K2as) -  (1+didi_1(2:K2as-1)) .* TguessS(2:K2as-1) + didi_1(2:K2as-1) .* TguessS(1:K2as-2)   ) ...
             ; ...
              -  Ttarget ...
                   .* one_over_lambda.^2 .*time ...
             ];       
    
%     checkJacob4 = max(abs(- Ttarget ...
%                    .* one_over_lambda.^2 .*tglobal))
%                
%    checkJacob5 = max(abs(- phi*dt/(hguess*hguess) ./(dx2(2:K2as-1) + dx2(1:K2as-2)) .* ...
%                     ( TguessS(3:K2as) -  (1+didi_1(2:K2as-1)) .* TguessS(2:K2as-1) + didi_1(2:K2as-1) .* TguessS(1:K2as-2)   ) ...
%               ))
               
    supdia2= [- phi*dt/(hguess*hguess) ./(dx2(2:K2as-1) + dx2(1:K2as-2)) .* ...
                    ( TguessS(3:K2as) -  TguessS(2:K2as-1)  )...
             ; ...
              sparse(K2-K2as,1)
             ];
                
    subdia2= [- phi*dt/(hguess*hguess) ./(dx2(2:K2as-1) + dx2(1:K2as-2)) .* ...
                    ( - didi_1(2:K2as-1) .* (TguessS(2:K2as-1) -  TguessS(1:K2as-2))   )...
             ; ...
              sparse(K2-K2as,1)
             ];     
               
  
    % multiply with partical derivatives of DT with respect to temperature
                
    diag   = diag   + diag2   .* dDTdTguessS(2:K2-1);
    supdia = supdia + supdia2 .* dDTdTguessS(3:K2);
    subdia = subdia + subdia2 .* dDTdTguessS(1:K2-2);
  
%     diag   = diag   + diag2   .* dDTdTguessS(2:K2as-1);
%     supdia = supdia + supdia2 .* dDTdTguessS(3:K2as);
%     subdia = subdia + subdia2 .* dDTdTguessS(1:K2as-2);
    
%     JTTsolINfull = myspdiag(subdia,diag,supdia);
% 
%     JTTsolIN = [JTTsolINfull,sparse(K2as-2,K2-K2as);...
%                 sparse(K2-K2as,K2as-1),speye(K2-K2as,K2-K2as),sparse(K2-K2as,1)];
%          
    JTTsolIN = myspdiag(subdia,diag,supdia);        
    % % ** ** **
%     % JTT 
%      
%     % Temperatur (solid) inner points: JTTsolIN (K2-2) x K2
%  
%     % (JTTsolid) Jacobian temperature residuals with respect to temperature - JTT in the solid phase (non-uniform grid) 
%     % maybe devide by ./xmesh(K1+2:K1+K2-1) (last row xmesh ->
%     % ones(K2-2,1);
%     
%     % NOTE: IN THE OLD IMPLEMENTATION OF THE DIAGONAL THE DIFFUSIVE TERM WAS PROB WRONG -
%     % not devided by xmesh !!!!!
%     
%     diag =  dx2(2:K2-1) ...
%             + phi*dt/(hguess*hguess) ./(dx2(2:K2-1) + dx2(1:K2-2)) .* ...
%                     ( DTguessS(3:K2) +  (1+didi_1(2:K2-1)) .* DTguessS(2:K2-1) + didi_1(2:K2-1) .* DTguessS(1:K2-2)   ) ...
%             + (1 - oldh/hguess) * phi * xmesh(K1+2:K1+K2-1);
% 
%     
% 
%     % superdiagonal
%     supdia= - phi*dt/(hguess*hguess) ./(dx2(2:K2-1) + dx2(1:K2-2)) .* ...
%                     ( DTguessS(3:K2) +  DTguessS(2:K2-1) ) ...
%             - (1 - oldh/hguess) * phi * xmesh(K1+2:K1+K2-1);
%  
%     
%     
%     %subdiagonal
%     subdia = - phi*dt/(hguess*hguess) ./(dx2(2:K2-1) + dx2(1:K2-2)) .* ...
%                     didi_1(2:K2-1) .* ( DTguessS(2:K2-1) +  DTguessS(1:K2-2) ); 
%     
% 
%        % ADD Coupling via the Diffusion constant
%     
%     
%     % Derivative with respect to DTguess
%     diag2  = - phi*dt/(hguess*hguess) ./(dx2(2:K2-1) + dx2(1:K2-2)) .* ...
%                     ( TguessS(3:K2) -  (1+didi_1(2:K2-1)) .* TguessS(2:K2-1) + didi_1(2:K2-1) .* TguessS(1:K2-2)   );
%     
%     supdia2= - phi*dt/(hguess*hguess) ./(dx2(2:K2-1) + dx2(1:K2-2)) .* ...
%                     ( TguessS(3:K2) -  TguessS(2:K2-1)  );
%                 
%     subdia2= - phi*dt/(hguess*hguess) ./(dx2(2:K2-1) + dx2(1:K2-2)) .* ...
%                     ( - didi_1(2:K2-1) .* (TguessS(2:K2-1) -  TguessS(1:K2-2))   );
%                
%   
%     % multiply with partical derivatives of DT with respect to temperature
%                 
%     diag   = diag   + diag2   .* dDTdTguessS(2:K2-1);
%     supdia = supdia + supdia2 .* dDTdTguessS(3:K2);
%     subdia = subdia + subdia2 .* dDTdTguessS(1:K2-2);
%   
%     
%     JTTsolIN = myspdiag(subdia,diag,supdia);

    % JTC vanishes no more due to coupling via the state dependent
    % diffusion constant - but: the concentration in the solid phase is
    % frozen in and not part of the state vector
    
    
    % THATS NOT TRUE - I skrewed it up.... there is ONE which does NOT
    % vanish There is one coupling of the TSolid(2) residual to the
    % DTsolid(1)
    
    JTCsolid = sparse(K2-2,1);
    
    % VANISHES NOW DUE TO IGNORED DEPENDENCE ON CONCENTRTION
%    JTCsolid(1,1) = subdia2(1) .* dDTdCguessS(1);

    % derivative of the 'Temp. Residuals'(liquid - inner points) with respect to h:
    % JThliquidIN (K1-2) x 1
    JThliqIN =   -(1-phi)/oldh              * ( ToldL(3:K1)   - ToldL(2:K1-1)   ) .* xmesh(2:K1-1) ...
                 - phi*oldh/(hguess*hguess) * ( TguessL(3:K1) - TguessL(2:K1-1) ) .* xmesh(2:K1-1)  ...
                 + 2.0*phi*dt/(hguess*hguess*hguess) ./(2.0*dx1) .* ...
                              (  TguessL(3:K1)   .* ( DTguessL(3:K1) + DTguessL(2:K1-1) ) ...
                               - TguessL(2:K1-1) .* ( DTguessL(3:K1) + 2.0 * DTguessL(2:K1-1) + DTguessL(1:K1-2) )  ...
                               + TguessL(1:K1-2) .* ( DTguessL(2:K1-1) + DTguessL(1:K1-2) ) ...
                              );
    

                          
% NEW

% (a) 2:K2as-1: original residual
% (b) K2-K2as: new 

    % derivative of the 'Temp. Residuals'(solid - inner points) with respect to h:
    % JThsolidIN (K2-2) x 1
    JThsolidINfull = -(1-phi)/oldh              * ( ToldS(3:K2as)   - ToldS(2:K2as-1)   ) .* xmesh(K1+2:K1+K2as-1) ...
                 - phi*oldh/(hguess*hguess) * ( TguessS(3:K2as) - TguessS(2:K2as-1) ) .* xmesh(K1+2:K1+K2as-1)  ...
                 + 2.0*phi*dt/(hguess*hguess*hguess) ./(dx2(2:K2as-1) + dx2(1:K2as-2)) .* ...
                              (  TguessS(3:K2as)   .* ( DTguessS(3:K2as) + DTguessS(2:K2as-1) ) ...
                               - TguessS(2:K2as-1) .* ( DTguessS(3:K2as) + DTguessS(2:K2as-1) + didi_1(2:K2as-1) .* ( DTguessS(2:K2as-1) + DTguessS(1:K2as-2) ) ) ...
                               + TguessS(1:K2as-2) .*                                       didi_1(2:K2as-1) .* ( DTguessS(2:K2as-1) + DTguessS(1:K2as-2) ) ...
                              );

% xeval = hguess .* xmesh(K1+K2as:K1+K2-1);
% 
% one_over_lambda = TintInvLam_global(xeval);
% 
% InitTemp = TintFct_global(xeval);
% 
% global Dconst_BBC
% Ttarget = InitTemp ...
%                   .* exp(Dconst_BBC .* one_over_lambda.^2 .*tglobal);
% % Ttarget = InitTemp ...
% %                   .* exp(DTguessS(K2as:K2-1) .* one_over_lambda.^2 .*tglobal);
% 
%      JThsolidINAs = xeval .* one_over_lambda .* Ttarget;
%JThsolidINAs = sparse(K2-K2as,1);
    
    dxtemp = min(diff(xeval));
    if(isempty(dxtemp))
       dxtemp = hguess; 
    end
    
    dlambda_dx = (   1./TintInvLam_global(xeval+0.5 * dxtemp)  ...
                   - 1./TintInvLam_global(xeval-0.5 * dxtemp) ...
                 ) ./dxtemp;
             
%     % CHECK
%      figure(54321)
% %     
%      plot(xeval,1./TintInvLam_global(xeval))
%      hold on
%      plot(xeval,dlambda_dx)
%      drawnow
%      pause
%              
    
    JThsolidINAs = xmesh(K1+K2as:K1+K2-1) .* Ttarget ...
        .*( one_over_lambda ...%;%...
          +  2.0 * DTguessS(K2as:K2-1)* time .* one_over_lambda.^3 .* dlambda_dx ...
          ) ;

%      checkJacob1 = max(abs( 2.0 * DTguessS(K2as:K2-1)* tglobal .* one_over_lambda.^3 .* dlambda_dx))
%      checkJacob2 = max(abs(JThsolidINAs)) 
%      checkJacob2a = max(abs(JThsolidINfull)) 
%      
%      size(JThsolidINfull)
%      size(JThsolidINAs)
%      K2
%      K2as
     JThsolidIN = [JThsolidINfull;JThsolidINAs];
                          
%     % derivative of the 'Temp. Residuals'(solid - inner points) with respect to h:
%     % JThsolidIN (K2-2) x 1
%     JThsolidIN = -(1-phi)/oldh              * ( ToldS(3:K2)   - ToldS(2:K2-1)   ) .* xmesh(K1+2:K1+K2-1) ...
%                  - phi*oldh/(hguess*hguess) * ( TguessS(3:K2) - TguessS(2:K2-1) ) .* xmesh(K1+2:K1+K2-1)  ...
%                  + 2.0*phi*dt/(hguess*hguess*hguess) ./(dx2(2:K2-1) + dx2(1:K2-2)) .* ...
%                               (  TguessS(3:K2)   .* ( DTguessS(3:K2) + DTguessS(2:K2-1) ) ...
%                                - TguessS(2:K2-1) .* ( DTguessS(3:K2) + DTguessS(2:K2-1) + didi_1(2:K2-1) .* ( DTguessS(2:K2-1) + DTguessS(1:K2-2) ) ) ...
%                                + TguessS(1:K2-2) .*                                       didi_1(2:K2-1) .* ( DTguessS(2:K2-1) + DTguessS(1:K2-2) ) ...
%                               );
%     
    


    
    % ** ** **
    % JCC 
    
    % Concentration (liquid) inner points: JCCliqIN (K1-2) x K1
 
    % (JCCliquid) Jacobian concentration residuals with respect to concentration - JCC in the liquid phase (uniform grid) 
    
        
    
    % diagonal
    diag =  dx1 ...
            + phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( DCguessL(3:K1) +  2.0 * DCguessL(2:K1-1) + DCguessL(1:K1-2)   ) ...
            + (1 - oldh/hguess) * phi * xmesh(2:K1-1);
        

    
    % superdiagonal
    supdia= - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( DCguessL(3:K1) +  DCguessL(2:K1-1) ) ...
            - (1 - oldh/hguess) * phi * xmesh(2:K1-1);

        
 
    
    
    % subdiagonal
    subdia = - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( DCguessL(2:K1-1) +  DCguessL(1:K1-2) ); 
    
    
  % ADD Coupling via the Diffusion constant
    
    
    % Derivative with respect to DCguess
    diag2  = - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( CguessL(3:K1) -  2.0 * CguessL(2:K1-1) + CguessL(1:K1-2)   );
    
    supdia2= - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( CguessL(3:K1) -  CguessL(2:K1-1) );
                
    subdia2= - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( -  CguessL(2:K1-1) +  CguessL(1:K1-2) ); 
             
    % multiply with partical derivatives of DC with respect to
    % concentration
    % NOTE: Dj with j the COLUMN index -> for the diag / sub / sup diag we
    % have to multiply with
    % diag2 = diag2 .* (2:K1-1) (i,j) 
    % supdia = supdia .* (3:K1) with j = i+1 (index of the diag
    % runsover i !!!!
    % subdia = subdia .* (1:K1-2) with j = i-1 (index of the diag
    % runsover i !!!!
    
    diag   = diag   + diag2   .* dDCdCguessL(2:K1-1);
    supdia = supdia + supdia2 .* dDCdCguessL(3:K1);
    subdia = subdia + subdia2 .* dDCdCguessL(1:K1-2);
      
    
    
    
    %JCCliqIN = spdiags([subdia diag supdia], 0:2, K1-2, K1);
    JCCliqIN = myspdiag(subdia, diag, supdia);


    % JCT vanishes no more due to coupling via the state dependent diffusion constant
    % dimension (K1-2) x K1 
    % Derivative with respect to DCguess multiplied with partical
    % derivatives of DC with respect to temperature
    
    diag   = diag2   .* dDCdTguessL(2:K1-1);
    supdia = supdia2 .* dDCdTguessL(3:K1);
    subdia = subdia2 .* dDCdTguessL(1:K1-2);
    
    
    JCTliqIN = myspdiag(subdia,diag,supdia);
    

    
    
    % JChliquidIN (K1-2) x 1
    % derivative of the 'Temp. Residuals'(liquid - inner points) with respect to h:


    JChliqIN =   -(1-phi)/oldh              * ( ColdL(3:K1)   - ColdL(2:K1-1)   ) .* xmesh(2:K1-1) ...
                 - phi*oldh/(hguess*hguess) * ( CguessL(3:K1) - CguessL(2:K1-1) ) .* xmesh(2:K1-1)  ...
                 + 2.0*phi*dt/(hguess*hguess*hguess) ./(2.0*dx1) .* ...
                              (  CguessL(3:K1)   .* ( DCguessL(3:K1) + DCguessL(2:K1-1) ) ...
                               - CguessL(2:K1-1) .* ( DCguessL(3:K1) + 2.0 * DCguessL(2:K1-1) + DCguessL(1:K1-2) )  ...
                               + CguessL(1:K1-2) .* ( DCguessL(2:K1-1) + DCguessL(1:K1-2) ) ...
                              );           
                              


   

        
else % Resolidification
% ** ** **
    % JTT 
    
    % Temperatur (liquid) inner points: JTTliqIN (K1-2) x K1
 
    % (JTTliquid) Jacobian temperature residuals with respect to temperature - JTT in the liquid phase (uniform grid) 

    
    % diagonal
    diag =  dx1 ...
            + phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( DTguessL(3:K1) +  2.0 * DTguessL(2:K1-1) + DTguessL(1:K1-2)   ) ...
            - (1 - oldh/hguess) * phi * xmesh(2:K1-1);


    
    % superdiagonal
    supdia= - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( DTguessL(3:K1) +  DTguessL(2:K1-1) );% ...
          %  - (1 - oldh/hguess) * phi * xmesh(2:K1-1);
 
    
    
    
    % subdiagonal
    subdia = - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( DTguessL(2:K1-1) +  DTguessL(1:K1-2) ) ... 
             + (1 - oldh/hguess) * phi * xmesh(2:K1-1);
 
    

    % ADD Coupling via the Diffusion constant
    
    
    % Derivative with respect to DTguess
    diag2  = - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( TguessL(3:K1) -  2.0 * TguessL(2:K1-1) + TguessL(1:K1-2)   );
    
    supdia2= - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( TguessL(3:K1) -  TguessL(2:K1-1) );
                
    subdia2= - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( -  TguessL(2:K1-1) +  TguessL(1:K1-2) ); 
             
    % multiply with partical derivatives of DT with respect to temperature
    % NOTE: Dj with j the COLUMN index -> for the diag / sub / sup diag we
    % have to multiply with
    % diag2 = diag2 .* (2:K1-1) (i,j) 
    % supdia = supdia .* (3:K1) with j = i+1 (index of the diag
    % runsover i !!!!
    % subdia = subdia .* (1:K1-2) with j = i-1 (index of the diag
    % runsover i !!!!
    
    diag   = diag   + diag2   .* dDTdTguessL(2:K1-1);
    supdia = supdia + supdia2 .* dDTdTguessL(3:K1);
    subdia = subdia + subdia2 .* dDTdTguessL(1:K1-2);
    

    %JTTliqIN = spdiags([subdia diag supdia], 0:2, K1-2, K1) ;
    JTTliqIN = myspdiag(subdia,diag,supdia);
    % first and last row are missing s. t. diagonal elements are at (i,i+1)


    % JTC vanishes no more due to coupling via the state dependent diffusion constant
    % dimension (K1-2) x K1 
    % Derivative with respect to DTguess multiplied with partical
    % derivatives of DT with respect to concentration
    
    diag   = diag2   .* dDTdCguessL(2:K1-1);
    supdia = supdia2 .* dDTdCguessL(3:K1);
    subdia = subdia2 .* dDTdCguessL(1:K1-2);
    
    
    JTCliqIN = myspdiag(subdia,diag,supdia);
    
    
    % ** ** **
    % JTT 
     
    % Temperatur (solid) inner points: JTTsolIN (K2-2) x K2
 
    % (JTTsolid) Jacobian temperature residuals with respect to temperature - JTT in the solid phase (non-uniform grid) 
    % maybe devide by ./xmesh(K1+2:K1+K2-1) (last row xmesh ->
    % ones(K2-2,1);
    
    % NOTE: IN THE OLD IMPLEMENTATION OF THE DIAGONAL THE DIFFUSIVE TERM WAS PROB WRONG -
    % not devided by xmesh !!!!!
    
    diag =  dx2(2:K2-1) ...
            + phi*dt/(hguess*hguess) ./(dx2(2:K2-1) + dx2(1:K2-2)) .* ...
                    ( DTguessS(3:K2) +  (1+didi_1(2:K2-1)) .* DTguessS(2:K2-1) + didi_1(2:K2-1) .* DTguessS(1:K2-2)   ) ...
            - (1 - oldh/hguess) * phi * didi_1(2:K2-1) .* xmesh(K1+2:K1+K2-1);

    

    % superdiagonal
    supdia= - phi*dt/(hguess*hguess) ./(dx2(2:K2-1) + dx2(1:K2-2)) .* ...
                    ( DTguessS(3:K2) +  DTguessS(2:K2-1) );% ...
           % - (1 - oldh/hguess) * phi * xmesh(K1+2:K1+K2-1);
 
    

    
    %subdiagonal
    subdia = - phi*dt/(hguess*hguess) ./(dx2(2:K2-1) + dx2(1:K2-2)) .* ...
                    didi_1(2:K2-1) .* ( DTguessS(2:K2-1) +  DTguessS(1:K2-2) ) ... 
             + (1 - oldh/hguess) * phi * didi_1(2:K2-1) .* xmesh(K1+2:K1+K2-1);
 
                

       % ADD Coupling via the Diffusion constant
    
    
    % Derivative with respect to DTguess
    diag2  = - phi*dt/(hguess*hguess) ./(dx2(2:K2-1) + dx2(1:K2-2)) .* ...
                    ( TguessS(3:K2) -  (1+didi_1(2:K2-1)) .* TguessS(2:K2-1) + didi_1(2:K2-1) .* TguessS(1:K2-2)   );
    
    supdia2= - phi*dt/(hguess*hguess) ./(dx2(2:K2-1) + dx2(1:K2-2)) .* ...
                    ( TguessS(3:K2) -  TguessS(2:K2-1)  );
                
    subdia2= - phi*dt/(hguess*hguess) ./(dx2(2:K2-1) + dx2(1:K2-2)) .* ...
                    ( - didi_1(2:K2-1) .* (TguessS(2:K2-1) -  TguessS(1:K2-2))   );
               
  
    % multiply with partical derivatives of DT with respect to temperature
                
    diag   = diag   + diag2   .* dDTdTguessS(2:K2-1);
    supdia = supdia + supdia2 .* dDTdTguessS(3:K2);
    subdia = subdia + subdia2 .* dDTdTguessS(1:K2-2);
  
    
    JTTsolIN = myspdiag(subdia,diag,supdia);

        % JTC vanishes no more due to coupling via the state dependent
    % diffusion constant - but: the concentration in the solid phase is
    % frozen in and not part of the state vector
    

    % THATS NOT TRUE - I skrewed it up.... there is ONE which does NOT
    % vanish There is one coupling of the TSolid(2) residual to the
    % DTsolid(1)
    
    JTCsolid = sparse(K2-2,1);
    % VANISHES NOW
%    JTCsolid(1,1) = subdia2(1) .* dDTdCguessS(1);


    
    % derivative of the 'Temp. Residuals'(liquid - inner points) with respect to h:
    % JThliquidIN (K1-2) x 1
    JThliqIN =   -(1-phi)/oldh              * ( ToldL(2:K1-1)   - ToldL(1:K1-2)   ) .* xmesh(2:K1-1) ...
                 - phi*oldh/(hguess*hguess) * ( TguessL(2:K1-1) - TguessL(1:K1-2) ) .* xmesh(2:K1-1)  ...
                 + 2.0*phi*dt/(hguess*hguess*hguess) ./(2.0*dx1) .* ...
                              (  TguessL(3:K1)   .* ( DTguessL(3:K1) + DTguessL(2:K1-1) ) ...
                               - TguessL(2:K1-1) .* ( DTguessL(3:K1) + 2.0 * DTguessL(2:K1-1) + DTguessL(1:K1-2) )  ...
                               + TguessL(1:K1-2) .* ( DTguessL(2:K1-1) + DTguessL(1:K1-2) ) ...
                              );
    
    

    % derivative of the 'Temp. Residuals'(solid - inner points) with respect to h:
    % JThsolidIN (K2-2) x 1
    JThsolidIN = -(1-phi)/oldh              * ( ToldS(2:K2-1)   - ToldS(1:K2-2)   ) .* didi_1(2:K2-1) .* xmesh(K1+2:K1+K2-1) ...
                 - phi*oldh/(hguess*hguess) * ( TguessS(2:K2-1) - TguessS(1:K2-2) ) .* didi_1(2:K2-1) .* xmesh(K1+2:K1+K2-1)  ...
                 + 2.0*phi*dt/(hguess*hguess*hguess) ./(dx2(2:K2-1) + dx2(1:K2-2)) .* ...
                              (  TguessS(3:K2)   .* ( DTguessS(3:K2) + DTguessS(2:K2-1) ) ...
                               - TguessS(2:K2-1) .* ( DTguessS(3:K2) + DTguessS(2:K2-1) + didi_1(2:K2-1) .* ( DTguessS(2:K2-1) + DTguessS(1:K2-2) ) ) ...
                               + TguessS(1:K2-2) .*                                       didi_1(2:K2-1) .* ( DTguessS(2:K2-1) + DTguessS(1:K2-2) ) ...
                              );


                              % ** ** **
    % JCC 
    
    % Concentration (liquid) inner points: JCCliqIN (K1-2) x K1
 
    % (JCCliquid) Jacobian concentration residuals with respect to concentration - JCC in the liquid phase (uniform grid) 
    
        
    % diagonal
    diag =  dx1 ...
            + phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( DCguessL(3:K1) +  2.0 * DCguessL(2:K1-1) + DCguessL(1:K1-2)   ) ...
            - (1 - oldh/hguess) * phi * xmesh(2:K1-1);
        

    
    % superdiagonal
    supdia= - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( DCguessL(3:K1) +  DCguessL(2:K1-1) );% ...
         %   - (1 - oldh/hguess) * phi * xmesh(2:K1-1);

        
 
    
    
    % subdiagonal
    subdia = - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( DCguessL(2:K1-1) +  DCguessL(1:K1-2) ) ... 
             + (1 - oldh/hguess) * phi * xmesh(2:K1-1);

    
  % ADD Coupling via the Diffusion constant
    
    
    % Derivative with respect to DCguess
    diag2  = - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( CguessL(3:K1) -  2.0 * CguessL(2:K1-1) + CguessL(1:K1-2)   );
    
    supdia2= - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( CguessL(3:K1) -  CguessL(2:K1-1) );
                
    subdia2= - phi*dt/(hguess*hguess) ./(2.0*dx1) .* ...
                    ( -  CguessL(2:K1-1) +  CguessL(1:K1-2) ); 
             
    % multiply with partical derivatives of DC with respect to
    % concentration
    % NOTE: Dj with j the COLUMN index -> for the diag / sub / sup diag we
    % have to multiply with
    % diag2 = diag2 .* (2:K1-1) (i,j) 
    % supdia = supdia .* (3:K1) with j = i+1 (index of the diag
    % runsover i !!!!
    % subdia = subdia .* (1:K1-2) with j = i-1 (index of the diag
    % runsover i !!!!
    
    diag   = diag   + diag2   .* dDCdCguessL(2:K1-1);
    supdia = supdia + supdia2 .* dDCdCguessL(3:K1);
    subdia = subdia + subdia2 .* dDCdCguessL(1:K1-2);
      
    
    
    
    %JCCliqIN = spdiags([subdia diag supdia], 0:2, K1-2, K1);
    JCCliqIN = myspdiag(subdia, diag, supdia);

    % JCT vanishes no more due to coupling via the state dependent diffusion constant
    % dimension (K1-2) x K1 
    % Derivative with respect to DCguess multiplied with partical
    % derivatives of DC with respect to temperature
    
    diag   = diag2   .* dDCdTguessL(2:K1-1);
    supdia = supdia2 .* dDCdTguessL(3:K1);
    subdia = subdia2 .* dDCdTguessL(1:K1-2);
    
    
    JCTliqIN = myspdiag(subdia,diag,supdia);

    % JChliquidIN (K1-2) x 1
    % derivative of the 'Temp. Residuals'(liquid - inner points) with respect to h:


    JChliqIN =   -(1-phi)/oldh              * ( ColdL(2:K1-1)   - ColdL(1:K1-2)   ) .* xmesh(2:K1-1) ...
                 - phi*oldh/(hguess*hguess) * ( CguessL(2:K1-1) - CguessL(1:K1-2) ) .* xmesh(2:K1-1)  ...
                 + 2.0*phi*dt/(hguess*hguess*hguess) ./(2.0*dx1) .* ...
                              (  CguessL(3:K1)   .* ( DCguessL(3:K1) + DCguessL(2:K1-1) ) ...
                               - CguessL(2:K1-1) .* ( DCguessL(3:K1) + 2.0 * DCguessL(2:K1-1) + DCguessL(1:K1-2) )  ...
                               + CguessL(1:K1-2) .* ( DCguessL(2:K1-1) + DCguessL(1:K1-2) ) ...
                              );           
                              


   

          
end

% Interface conditions at the moving fron
% (a) energy conservation
JBCmovA = sparse(1,Z+K1+2) ;


prefac = dt*Aparam*phi/hguess;
% Energy flux in liquid
JBCmovA(1,K1)  =  prefac * ( 1.5/dx1 * KappaThermguessL ...
                  + (1.5 * TguessL(K1) - 2.0     * TguessL(K1-1) + 0.5 * TguessL(K1-2) )/dx1 * dKappaThermdTguessL );
JBCmovA(1,K1-1)=  prefac * (-2.0/dx1 * KappaThermguessL );
JBCmovA(1,K1-2)=  prefac * (0.5/dx1 * KappaThermguessL);
% Energy flux in solid
JBCmovA(1,K1+1)=  prefac * (     dwa   * KappaThermguessS ...
                  + (dwb * TguessS(3) - (dwa+dwb)* TguessS(2)    + dwa *   TguessS(1)  )     * dKappaThermdTguessS );
JBCmovA(1,K1+2)=  prefac * (-(dwa+dwb) * KappaThermguessS );
JBCmovA(1,K1+3)=  prefac * (   dwb     * KappaThermguessS);

% additional coupling to the concentration at the interface (via the solute
% con. dependence of thermal diffusivity)
JBCmovA(1,Z+K1)  =  prefac * (1.5 * TguessL(K1) - 2.0     * TguessL(K1-1) + 0.5 * TguessL(K1-2) )/dx1 * dKappaThermdCguessL;
JBCmovA(1,Z+K1+1)=  prefac * (dwb * TguessS(3) - (dwa+dwb)* TguessS(2)    + dwa *   TguessS(1)  )     * dKappaThermdCguessS;

% derivitave with respect to h - index Z+K1+2
JBCmovA(1,Z+K1+2)= 1 - prefac/hguess *( (1.5 * TguessL(K1) - 2.0     * TguessL(K1-1) + 0.5 * TguessL(K1-2) )/dx1 * KappaThermguessL ...
                                      - (-dwb * TguessS(3) +(dwa+dwb)* TguessS(2)    - dwa *   TguessS(1)  )     * KappaThermguessS  );

                                  
   
                                  
                                  
% (b) solute conservation 
JBCmovB = sparse(1,Z+K1+2) ;


JBCmovB(1,Z+K1) =   1.5*phi*dt*DSguess + ...
                    dt * phi * ( 1.5 * CguessL(K1) - 2.0 * CguessL(K1-1) + 0.5 * CguessL(K1-2) ) * dDCdCguessL(K1) ... 
                  + 0.5*phi*dx1 * (hguess*hguess - oldh*oldh);

JBCmovB(1,Z+K1-1) = -2.0*phi*dt*DSguess;
JBCmovB(1,Z+K1-2) =  0.5*phi*dt*DSguess;
JBCmovB(1,Z+K1+1) = -0.5*phi*dx1 * (hguess*hguess - oldh*oldh);

% Temperature via T dependence of the diffusion constant
JBCmovB(1,K1)     = dt * phi * ( 1.5 * CguessL(K1) - 2.0 * CguessL(K1-1) + 0.5 * CguessL(K1-2) ) * dDCdTguessL(K1);

% h derivative
JBCmovB(1,Z+K1+2) =  hguess*dx1* ...
                    (  phi    *  (CguessL(K1) - CguessS(1) )  ...
                    + (1-phi) *  (ColdL(K1)   - ColdS(1)   )  );


                
%(c) continuity of T - index K1+1
JBCmovC = sparse(1,Z+K1+2) ;
JBCmovC(1,K1) = 1;
JBCmovC(1,K1+1) = -1;
                
% INTERFACE RESPONSE FUNCTIONS
% (a) Temperature as a fct. of front velocity - in the documentation
% \tilde{f}

                
%JIRFmovT = sparse(1,Z+K1+2);       

% function of hdot, interface Temp, concentrations a interface and the
% state (amorphous / crystaline) of the solid
JIRFmovT = getJacobIRFmovT(oldh,hguess,dt,phi,ToldL(K1),TguessL(K1),ColdL(K1),CguessL(K1),ColdS(1),CguessS(1), stateSold,stateSguess);       

 
if(~hincrease)               
% (b) jump in the concentration - in the documentation
% \tilde{g}  (only relevant for resolidification part - \dot{h} < 0
% function of hdot, interface Temp, concentrations a interface and the
% state (amorphous / crystaline) of the solid

%JIRFmovC = sparse(1,Z+K1+2);
JIRFmovC = getJacobIRFmovC(oldh,hguess,dt,phi,ToldL(K1),TguessL(K1),ColdL(K1),CguessL(K1),ColdS(1),CguessS(1), stateSold,stateSguess);


else
% (b2) for \dot{h} >0 (melting) we instead keep the solid concentration
% value at the moving interface at whatever value of clab(h) -
% interpolation

JICmovC = sparse(1,Z+K1+2);
JICmovC(1,Z+K1+1) = 1;

% NEW Derivative with respect to hguess
  dxtemp = 1e-3*hguess;
  dCSenfore_dh = (  getInterfaceCSfromCSlab(hguess +0.5 * dxtemp) ...
                  - getInterfaceCSfromCSlab(hguess -0.5 * dxtemp) ...
                 ) ./dxtemp; 
  
             
  JICmovC(1,Z+K1+2) = - dCSenfore_dh;

% ATTENTION - there is an implicit dependencnce on h due to the changing
% location where clab is interpolated at. Thus is not analytically
% available so we have to either do finite differences: - d(cenforce)/dh \aprox -
% (CSenforceGuess - CSenforceOld)/(hguess-oldh)
%JICmovC(1,Z+K1+2) = -(CSenforceGuess - CSenforceOld)./(hguess-oldh+1e-14);

% OR: because that does not really work due to problems with the
% interpolation.... we keep the enforced concentration fixed - see
% concSolidEnforceGlobal which is updated in the main loop and OUTSIDE the
% Newton iteration
end

% Boundary conditions at surface (x=0)
% (a) Temperature at the surface - in the documentation \tilde{F} (capital
% F)
% depends on TguessL(2), TguessL(1), hguess, tglobal (explicit time
% dependence 

%JBCsurfT = sparse(1,Z+K1+2);
JBCsurfT = getJacobBCsurfT(TguessL(1),TguessL(2),hguess,tglobal);



% (b) Concentration at the surface - in the documentation \tilde{G} (capital G)
% NOTE: Can be explicitly time dependent - provide time as a parameter -
% not as an element of the state vector
% depends on CguessL(2), CguessL(1), Tguess(1) - surface temperature,
% hguess, tglobal

%JBCsurfC = sparse(1,Z+K1+2);
JBCsurfC = getJacobBCsurfC(CguessL(1),CguessL(2),TguessL(1),hguess,tglobal);




% Boundary conditions in the bulk x \to \infty (actually L2)
% Temperature approaches bulk temperature (=0 in our units). Alternatively
% we can consider an effective BC to account fo the finite value of L2 - in
% the documentation \tilde{H}

%JBCbulkT = sparse(1,Z+K1+2);
%JBCbulkT(1,Z) = 1.0;

if (hincrease)
JBCbulkT = getJacobBCbulkT(TguessS(K2),hguess,tglobal);
else
JBCbulkT = getJacobBCbulkTresold(TguessS(K2),hguess,tglobal,ToldS(K2-3));
end

% Put together the jacobian


if (hincrease)

jacob = [JBCsurfT; ...
        JTTliqIN,sparse(K1-2,K2),JTCliqIN,sparse(K1-2,1),JThliqIN; ...
        JBCmovA; ...
        JBCmovC; ...
%        sparse(K2-2,K1),JTTsolIN, sparse(K2-2,K1+1),JThsolidIN;...
        sparse(K2-2,K1),JTTsolIN, sparse(K2-2,K1),JTCsolid,JThsolidIN;...
        JBCbulkT; ...
        JBCsurfC; ...
        JCTliqIN,sparse(K1-2,K2),JCCliqIN,sparse(K1-2,1),JChliqIN; ...
        JBCmovB; ...
        JICmovC; ...
        JIRFmovT];
else
jacob = [JBCsurfT; ...
        JTTliqIN,sparse(K1-2,K2),JTCliqIN,sparse(K1-2,1),JThliqIN; ...
        JBCmovA; ...
        JBCmovC; ...
%        sparse(K2-2,K1),JTTsolIN, sparse(K2-2,K1+1),JThsolidIN;...
        sparse(K2-2,K1),JTTsolIN, sparse(K2-2,K1),JTCsolid,JThsolidIN;...
        JBCbulkT; ...
        JBCsurfC; ...
        JCTliqIN,sparse(K1-2,K2),JCCliqIN,sparse(K1-2,1),JChliqIN; ...
        JBCmovB; ...
%        JICmovC; ...
        JIRFmovC;
        JIRFmovT];
    
    
end
   

% CHECK

% jacobNum = getjacobNum(uold,uguess,dt);
% 
% size(jacob)
% size(jacobNum)
% checkJacob = jacob - jacobNum;
% %figure(432524)
% %imagesc((abs(checkJacob)));
% 
% figure(432523)
% imagesc(log(abs(checkJacob)));
% % figure(432522)
% % imagesc(log(abs(jacob(:,end))));
% % figure(432521)
% % imagesc(log(abs(jacobNum(:,end))));
% 
% checkerror = max(max(abs(checkJacob)))
% 
% [maxA,ind] = max(abs(checkJacob(:)));
% [m,n] = ind2sub(size(checkJacob),ind)
% 
% checkerrorrel=max(max(abs(checkJacob./jacobNum)))
% jacob(K1+K2as:K1+K2as+4,end)
% jacobNum(K1+K2as:K1+K2as+4,end)
% checkJacob(K1+K2as-1:K1+K2,end)
% size(checkJacob(K1+K2as-1:K1+K2,end))

% if(hguess*xmesh(end) >12000)
% %   pause
%    
%    
%    
%    jacobNum = getjacobNum(uold,uguess,dt);
% % 
% % size(jacob)
% % size(jacobNum)
%  checkJacob = jacob - jacobNum;
%  %figure(432524)
%  %imagesc((abs(checkJacob)));
% % 
%  figure(432523)
%  imagesc(log(abs(checkJacob)));
% % % figure(432522)
% % % imagesc(log(abs(jacob(:,end))));
% % % figure(432521)
% % % imagesc(log(abs(jacobNum(:,end))));
% % 
%  checkerror = max(max(abs(checkJacob)))
% % 
%  [maxA,ind] = max(abs(checkJacob(:)));
%  [m,n] = ind2sub(size(checkJacob),ind)
% % 
%  checkerrorrel=max(max(abs(checkJacob./jacobNum)))
% % jacob(K1+K2as:K1+K2as+4,end)
% % jacobNum(K1+K2as:K1+K2as+4,end)
% % checkJacob(K1+K2as-1:K1+K2,end)
% % size(checkJacob(K1+K2as-1:K1+K2,end))
%  
% %pause
% end

%pause
        
end
