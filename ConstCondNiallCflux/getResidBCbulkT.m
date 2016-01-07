function [ResBCbulkT] = getResidBCbulkT(Tguess,hguess,tglobal)


global Dconst_BBC thresh_check_faraway_BBC TempThreshold_BBC
global L2
global TintFct_global TintInvLam_global


% Check if we are far enough away from the surface

% Three checks for validity of the asymptotics (see notes 08/12/12)
% (1) lambda >= sigma approx 2
% (2) short time: t < lambd^2/D
% (3) for field: L2h(t) >> sqrt(4 D t)

% (this is for check (3)
check=hguess./sqrt(tglobal);

if (check < thresh_check_faraway_BBC)
   
    check
    error('residual: Asymptotics for bulk T: evaluating too close to surface')
    
end

xeval = hguess*L2;
%InitTemp = getInitTempFct(xeval,Tsurf_BBC, Ldec_BBC);
InitTemp = TintFct_global(xeval);

Tthreshold = TempThreshold_BBC;

    if(InitTemp < Tthreshold)
         Ttarget = 0.0;
    else
    % NOTE: This will become lambda dependent
    %    if(tglobal < tasymp_BBC)
%         Ttarget = InitTemp ...
%                   * exp(Dconst_lambdasq_BBC * tglobal);
%        Ttarget = InitTemp ...
%                  * exp(Dconstpar_gaussian_BBC *hguess*hguess * tglobal);
    
    % CHECKS
    
    % (1)
      %   one_over_lambda = getInitTempLambdaFct(xeval);
       one_over_lambda = TintInvLam_global(xeval);
    
        if (1./one_over_lambda < 2)
         %    one_over_lambda
         %    InitTemp
           % error('residual: Asymptotics for bulk T: lambda too small')
           fprintf('warning: residual: Asymptotics for bulk T: lambda very small')
           lambda = 1./  one_over_lambda
           
   
        end
    %  (2)  
        t_threshold = 1./ one_over_lambda ./one_over_lambda./Dconst_BBC;

        if (tglobal > t_threshold)
            tglobal
            t_threshold
            error('residual: Asymptotics for bulk T: time too late')
        end
            
        Ttarget = InitTemp ...
                  * exp(Dconst_BBC .* one_over_lambda .* one_over_lambda .*tglobal);
    
     end
%Ttarget
ResBCbulkT = Tguess-Ttarget;

    
return