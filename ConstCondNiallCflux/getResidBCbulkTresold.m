function [ResBCbulkT] = getResidBCbulkTresold(dt,oldh,hguess,TguessS,t_global,ToldCheck)

global Tau_BBC deltaTau_BBC T0fit_BBC x0_BBC lambdafit_BBC Dconst_BBC xinftyTau_deltaTau_BBC curveparam_BBC TempThreshold_BBC plotflag_BBC

global K1 K2 Z L2 xmesh


beta_param = 0.5;%   use the interpolation in the intervall [tau, tau + beta * deltaTau]
% new structure:


TempThreshold = TempThreshold_BBC;%1e-12;
if (ToldCheck < TempThreshold)
    ResBCbulkT = TguessS(K2);
    return
else    
xeval = hguess * L2;

if (Tau_BBC <= t_global && t_global <Tau_BBC + beta_param * deltaTau_BBC && xeval > xinftyTau_deltaTau_BBC)
    
%    Ttarget = T0fit_BBC .* exp( - (xeval - x0_BBC)./lambdafit_BBC ) ...
%              .* exp(Dconst_BBC ./ lambdafit_BBC ./lambdafit_BBC  .*(t_global - Tau_BBC))

    
    TimeDependence = exp(Dconst_BBC ./ lambdafit_BBC ./lambdafit_BBC  .*(t_global - Tau_BBC));
    
 %   Ttau = T0fit_BBC .* exp( - (xeval - x0_BBC)./lambdafit_BBC )
    Ttau = T0fit_BBC .* exp( - (xeval - x0_BBC)./lambdafit_BBC + curveparam_BBC.*(xeval - x0_BBC).^2 );
    Ttarget = Ttau .* TimeDependence;
    
          
 %   figure(129)
    %subplot(1,2,1)
 %   hold on
 %   plot(xeval - x0_BBC,log(Ttarget),'xr')
   
%    figure(130)
%    subplot(1,2,2)
%hold on
%    plot(t_global-Tau_BBC,TimeDependence-1.0,'x')
    
    
    
    
          
else
    
   % renew the approximation
    Tau_BBC = t_global;

    %deltaTau_BBC = 10;
    
    deltaX = 200;

    
    xinftyTau = L2.* hguess;
    
    
    for counter =1:10
    xinftyTau_deltaTau_BBC = xinftyTau - deltaX;
    
   
    % find all values of the grid in the intervall at time tau!!!!
    xinfty_later_index = find(hguess*xmesh(K1+1:Z)> xinftyTau_deltaTau_BBC, 1,'first');

    
    % now we fit an exponential to the part of the temperature field: 



    % fit temperature in this range to an exponential
    % TRY THAT
    %x0_BBC = xinftyTau_deltaTau_BBC;
    x0_BBC = xinftyTau;
    p = polyfit(hguess*xmesh(K1+xinfty_later_index:K1+K2-5) - x0_BBC,log(TguessS(xinfty_later_index:K2-5)),2);

     
%     lambdafit_BBC = -1./p(1);
%     T0fit_BBC    = exp(p(2)); % that's 
    curveparam_BBC  = p(1);
    lambdafit_BBC = -1./p(2);
    T0fit_BBC    = exp(p(3)); % that's 
    curvature_lengthscale = sqrt(1./abs(p(1)));
   
    %pause
    % check if detlaTau is small enough
    % (1) validity of the exponential fit
    % (2) time short enough
    paramlarge = 10;
    parammoderate = 2;
    
    deltaTau_BBC = 1./paramlarge * lambdafit_BBC * lambdafit_BBC ./Dconst_BBC ;

    
    check1 = K2 - xinfty_later_index > 8; % enough points for a fit
%     check1e = K2 - xinfty_later_index > 4 % enough points for a fit
    
    % if there are two few points the fit is arbitrary so that the rest of
    % the checks does not 'work'
%     if (~check1e)
%          deltaX = 2*deltaX;
%          continue
%     end
%     
    
    check2a = deltaX < curvature_lengthscale;
    %check2a = deltaX <1./parammoderate * curvature_lengthscale;
    %check2a = deltaX <1./paramlarge * curvature_lengthscale

    check2b = K2 - xinfty_later_index < 35; % too many points for a fit
    
    % too large
    check2 = check2a && check2b;
    
%     if(ToldCheck < parammoderate.*TempThreshold)
%     check2 = true
%     end
%     
    
    
    
   
    if (check1 && check2)
   
         if (plotflag_BBC)
            figure(129)
            hold off
            plot(hguess*xmesh(K1+xinfty_later_index:K1+K2-3) - x0_BBC,log(TguessS(xinfty_later_index:K2-3)),'o')
            hold on
            tempx=hguess*xmesh(K1+xinfty_later_index:K1+K2-3) - x0_BBC;
            plot(tempx, p(3) + p(2)*tempx + p(1) * tempx.*tempx,'k')
            plot(tempx, p(3) + p(2)*tempx,'r')
            title('Resolidification: Exponential tails - fit with (x0,T0) = (L2*h,T(h))');
            xlabel('x - x0');
            ylabel('log(T - T0)');
   
         end
   
        
        
        break
        
    elseif (check1 && ~check2)
        
        deltaX = 0.5*deltaX;
    
    %elseif (~check1 && check2)
    elseif (~check1)
    
        deltaX = 2*deltaX;
        
        
    end
   
    
    end

    if(~check1 && ~check2)
        deltaX
        curvature_lengthscale
         if (plotflag_BBC)
            figure(129)
            hold off
            plot(hguess*xmesh(K1+xinfty_later_index:K1+K2-3) - x0_BBC,log(TguessS(xinfty_later_index:K2-3)),'o')
            hold on
            tempx=hguess*xmesh(K1+xinfty_later_index:K1+K2-3) - x0_BBC;
            plot(tempx, p(3) + p(2)*tempx + p(1) * tempx.*tempx,'k')
            plot(tempx, p(3) + p(2)*tempx,'r')
            title('Resolidification: Exponential tails - fit with (x0,T0) = (L2*h,T(h))');
            xlabel('x - x0');
            ylabel('log(T - T0)');
   
        end
        
        error('Can not find a suitable deltaX')
    end
    
    
    check4 = lambdafit_BBC < parammoderate.* curvature_lengthscale;
    % only then the expoenential approximation hold
     if (~check4)
       %error('check4 failed') 
       fprintf('check4 failed') 
   end

   %pause
        
       TimeDependence = exp(Dconst_BBC ./ lambdafit_BBC ./lambdafit_BBC  .*(t_global - Tau_BBC));
    
 %   Ttau = T0fit_BBC .* exp( - (xeval - x0_BBC)./lambdafit_BBC )
    Ttau = T0fit_BBC .* exp( - (xeval - x0_BBC)./lambdafit_BBC + curveparam_BBC.*(xeval - x0_BBC).^2 );
    Ttarget = Ttau .* TimeDependence;
       
    
    
end

ResBCbulkT = TguessS(K2)-Ttarget;


end
return