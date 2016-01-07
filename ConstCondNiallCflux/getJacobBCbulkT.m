function [JBCbulkT] = getJacobBCbulkT(Tguess,hguess,tglobal);

global K1  L2 Z 
global Dconst_BBC  TempThreshold_BBC
global TintFct_global TintInvLam_global
% Vector: T[1:K1] liquid, T[K1+1,K1+K2] solid, C[1:K1] liquid, C[K1+1] (SINGLE
% POINT) solid, h

% Boundary conditions at surface (x=0)
% (a) Temperature at the surface - in the documentation \tilde{F} (capital
% F)
% depends on TguessL(2), TguessL(1), hguess, tglobal (explicit time
% dependence 

JBCbulkT = sparse(1,Z+K1+2);

% Derivative w.r.t Tguess

JBCbulkT(1,Z) = 1.0;

% Derv w.r.t. hguess
xeval = hguess*L2;
%InitTemp = getInitTempFct(xeval,Tsurf_BBC, Ldec_BBC);
InitTemp = TintFct_global(xeval);
Tthreshold = TempThreshold_BBC;%1e-12;

if(InitTemp < Tthreshold)
    %Ttarget = 0.0;
    % nothing - derivatives vanishes
else    
%    if(tglobal < tasymp_BBC)
%         Ttarget = InitTemp ...
%                   * exp(Dconst_lambdasq_BBC * tglobal);
       % Ttarget = InitTemp ...
       %           * exp(Dconstpar_gaussian_BBC *hguess*hguess * tglobal);
    
       
       
    %one_over_lambda = getInitTempLambdaFct(xeval);
    one_over_lambda = TintInvLam_global(xeval);
    
    dxtemp = 1e-4*xeval;
    
    dlambda_dx = (   1./TintInvLam_global(xeval+0.5 * dxtemp)  ...
                   - 1./TintInvLam_global(xeval-0.5 * dxtemp) ...
                 ) ./dxtemp;
    
     Ttarget = InitTemp ...
                  * exp(Dconst_BBC .* one_over_lambda .* one_over_lambda .*tglobal);
    
    %     error     
    %JBCbulkT(1, Z+K1+2) = + 2.* L2 * xeval ./Ldec_BBC ./Ldec_BBC .*Ttarget;  
    %JBCbulkT(1, Z+K1+2) = L2 .* one_over_lambda .*Ttarget  ;
    
    JBCbulkT(1, Z+K1+2) = L2 .* Ttarget ...
        .*( one_over_lambda ...%;%...
          +  2.0 * Dconst_BBC * tglobal .* one_over_lambda.^3 .* dlambda_dx ...
          ) ;
    
%    normJacob = abs(L2.* one_over_lambda .*Ttarget) 
    %error
%    else
%        InitTemp
%        error('jacobian; Asymptotics for bulk T: time to large for approx. to be valid and initial Temp. too large ')
%    
%    end
end
%JBCbulkT
return

