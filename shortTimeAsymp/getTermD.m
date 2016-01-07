function [output] = getTermD(xeval, U_tilde, t_tilde, DiffConstant,alphaParam)
global plotflagIC quadgkRelTol_STA
    
    f = @(tauint) getIntFctD1(tauint,xeval,U_tilde, t_tilde, DiffConstant);
    sqr_t_tilde = sqrt(t_tilde);
    if plotflagIC == true
         figure(2);
        subplot(2,1,2)
        %figure(4)
        plot(linspace(0,sqr_t_tilde,100),f(linspace(0,sqr_t_tilde,100)));
        title('Integrands Term D');
        xlabel('time');
    end
    
    %int1 = quad(f,0,sqr_t_tilde);
  %  int1 = quadgk(f,0,sqr_t_tilde);
    int1 = quadgk(f,0,sqr_t_tilde,'RelTol',quadgkRelTol_STA);
    
    g = @(tauint) getIntFctD2(tauint,xeval,U_tilde, t_tilde, DiffConstant);
   
    if plotflagIC == true
        hold on;
        plot(linspace(0,sqr_t_tilde,100),g(linspace(0,sqr_t_tilde,100)));
        hold off;
    end
    
   % int2 = quad(g,0,sqr_t_tilde);
   % int2 = quadgk(g,0,sqr_t_tilde);
    int2 = quadgk(g,0,sqr_t_tilde,'RelTol',quadgkRelTol_STA);
    

output = -alphaParam .* U_tilde ./ sqrt(DiffConstant*pi) ...
         .* exp(-U_tilde .* U_tilde .* t_tilde .* t_tilde .* 0.25./ DiffConstant) ...
         .* ( ...
              exp(0.5 ./DiffConstant .* (U_tilde .* t_tilde - xeval) .* U_tilde) ...
                  .* int1 ...
            + exp(0.5 ./DiffConstant .* (U_tilde .* t_tilde + xeval) .* U_tilde) ...
                  .* int2 ...
            );      
     
    

return

function [integrand] = getIntFctD1(tau,x,U,ttild,D)


    integrand =     exp( - 0.25./D .* (U.*ttild - x) .* (U.*ttild - x) ...
                                   ./ tau ./tau );
                               
return


function [integrand2] = getIntFctD2(tau,x,U,ttild,D)


    integrand2 =     exp( - 0.25./D .* (U.*ttild + x) .* (U.*ttild + x) ...
                                   ./ tau ./tau );
                               
return
