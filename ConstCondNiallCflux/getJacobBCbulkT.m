function [JBCbulkT] = getJacobBCbulkT(Tguess,hguess,tglobal);

global K1 L2 Z 
global  TempThreshold_BBC
global TintFct_global TintInvLam_global
% Vector: T[1:K1] liquid, T[K1+1,K1+K2] solid, C[1:K1] liquid, C[K1+1] (SINGLE
% POINT) solid, h

% Boundary conditions far into bulk during melting


JBCbulkT = sparse(1,Z+K1+2); %intialize 

DT_BBC = getDTS(Tguess,2);
DT_BBCdT = getDTSdT(Tguess,2);

% calculate position we are evaluating at and initial temperature
% parameters
xeval = hguess*L2;
InitTemp = TintFct_global(xeval);
Tthreshold = TempThreshold_BBC;%1e-12;


one_over_lambda = TintInvLam_global(xeval);
Ttarget = InitTemp ...
        * exp(DT_BBC .* one_over_lambda .* one_over_lambda .*tglobal);
    

% Derivative w.r.t Tguess
JBCbulkT(1,Z) = (1.0 - tglobal*Ttarget*DT_BBCdT*one_over_lambda*one_over_lambda);

% Derv w.r.t. hguess
if(InitTemp < Tthreshold)

    % nothing - derivatives vanishes
else    

    dxtemp = 1e-4*xeval;
    
    dlambda_dx = (   1./TintInvLam_global(xeval+0.5 * dxtemp)  ...
                   - 1./TintInvLam_global(xeval-0.5 * dxtemp) ...
                 ) ./dxtemp;

    JBCbulkT(1, Z+K1+2) = L2 .* Ttarget ...
        .*( one_over_lambda ...%;%...
          +  2.0 * DT_BBC * tglobal .* one_over_lambda.^3 .* dlambda_dx ...
          ) ;

end

return

