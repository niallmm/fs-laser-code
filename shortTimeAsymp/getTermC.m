function [output] = getTermC(x,t_tilde,DiffConstant,Linfty,maxRelError)
    global plotflagIC
global TintFct_global quadgkRelTol_STA
% Integrate 
% f = @(xint) getInitTempFct(xint, Tsurface,Ldecay) .* ...
%                 getGreensFct(x, xint, t_tilde, 0.0,DiffConstant);
   
    f = @(xint) TintFct_global(xint) .* ...
                getGreensFct(x, xint, t_tilde, 0.0,DiffConstant);
    
   % ezplot(@(x) f(x),0,Linfty,fig);
    if plotflagIC == true
         figure(2);
        subplot(2,1,1)
       %figure(3);
         plot(linspace(0,Linfty,100),f(linspace(0,Linfty,100)));
        title('Integrand Term C');
        xlabel('space');
    end
   % output = quad(f,0,Linfty);
%    output = quadgk(f,0,Linfty);
    output = quadgk(f,0,Linfty,'RelTol',quadgkRelTol_STA);
    
   % interror = quad(f,Linfty,2.0*Linfty) ./ output;
    interror = quadgk(f,Linfty,2.0*Linfty) ./ output;

    if abs(interror) > maxRelError
        'Integration error too large: ',interror
        error('increase Linfty') 
    end

            

end
