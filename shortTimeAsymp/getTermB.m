function [output] = getTermB(U_tilde,t_tilde,DiffConstant,alphaParam)
    global plotflagIC quadgkRelTol_STA
    % Integrate 
    %f = @(tint) getGreensFct(U_tilde.*t_tilde, U_tilde.*tint, t_tilde, tint,DiffConstant);
    
    %t_tilde
    %f(0.55)
    %f(0.59)
    %%f(0.59999)
    %f(t_tilde (1.0 - ))
    %ezplot(@(t) f(t),0.1,t_tilde);
    %plot(linspace(0,t_tilde),f(linspace(0,t_tilde)));
    
    %output = quad(f,0,t_tilde-realmin);

    % PROBLEM: Integrable Singularity at tint = t_tilde (due to
    % 1/sqrt(t_tile-tint) term
    
    % Solution: change variable to remove singularity
    

    f = @(tauint) getIntFctB(tauint,U_tilde, t_tilde, DiffConstant);
    sqr_t_tilde = sqrt(t_tilde);

    if plotflagIC == true
        figure(1);
        subplot(2,1,2)
        %title('Integrand Term B');
        %figure(2)
        plot(linspace(0,sqr_t_tilde,100),f(linspace(0,sqr_t_tilde,100)));
        title('Integrand Term B');
        xlabel('time');
    end
    %output = -alphaParam .* U_tilde ./ sqrt(DiffConstant*pi) ...
    %        .* quad(f,0,sqr_t_tilde);
    output = -alphaParam .* U_tilde ./ sqrt(DiffConstant*pi) ...
            .* quadgk(f,0,sqr_t_tilde,'RelTol',quadgkRelTol_STA);
    
end

function [integrand] = getIntFctB(tau,U,ttild,D)

    % NOTE: Make that faster - vecortized treatment of tau = 0
    Usq_D = U.*U./D;
    
    expTerm = exp(Usq_D .*ttild);
    
%    if  tau == 0
        
%    else
    integrand =   exp( -0.25 * Usq_D .* tau .* tau  ) ...
                    .* (   1.0   ...
                         + expTerm .* exp(-Usq_D .* ttild .*ttild ./tau ./tau ) ...
                       );
%    end
    


end
