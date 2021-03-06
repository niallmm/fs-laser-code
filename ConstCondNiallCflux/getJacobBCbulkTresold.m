function [JBCbulkT] = getJacobBCbulkTresold(Tguess,hguess,t_global, ToldCheck);

global K1  L2 Z 
%global Tsurf_BBC Ldec_BBC tasymp_BBC Dconstpar_gaussian_BBC
global Tau_BBC deltaTau_BBC T0fit_BBC x0_BBC lambdafit_BBC Dconst_BBC xinftyTau_deltaTau_BBC curveparam_BBC  TempThreshold_BBC
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
TempThreshold = TempThreshold_BBC;%1e-12;

if (ToldCheck < TempThreshold)
    %Ttarget = 0.0;
    % nothing - derivatives vanishes
else
xeval = hguess*L2;

% DON'T check the approximation validity - that is done in the residulal
% and then the approx may be updated
beta_param = 0.5;%   use the interpolation in the intervall [tau, tau + beta * deltaTau]

    if (Tau_BBC <= t_global && t_global <Tau_BBC + beta_param * deltaTau_BBC && xeval > xinftyTau_deltaTau_BBC)
   
 
       TimeDependence = exp(Dconst_BBC ./ lambdafit_BBC ./lambdafit_BBC  .*(t_global - Tau_BBC));
    
 %   Ttau = T0fit_BBC .* exp( - (xeval - x0_BBC)./lambdafit_BBC )
    Ttau = T0fit_BBC .* exp( - (xeval - x0_BBC)./lambdafit_BBC + curveparam_BBC.*(xeval - x0_BBC).^2 );
    Ttarget = Ttau .* TimeDependence;
       

             
    %JBCbulkT(1, Z+K1+2) = + 2.* L2 * xeval ./Ldec_BBC ./Ldec_BBC .*Ttarget;  
    JBCbulkT(1, Z+K1+2) = L2 ./lambdafit_BBC .*Ttarget  ;
    
    %normJacob = abs(L2 ./lambdafit_BBC .*Ttarget)
    %error
    
    else
        xeval
        xinftyTau_deltaTau_BBC
       error('getJacobBCbulkTresid: outside valid region') 
    end
    
end
return

