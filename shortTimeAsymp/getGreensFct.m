function [output] = getGreensFct(x,x0,t,t0,D)
    % Diffusion constant
    %D = 1.0;
    
    t_t0 = t - t0;
    inv4D_t_t0 = 1./(4*D*t_t0); 
    
    output  = (     exp( - (x - x0).*(x - x0).*inv4D_t_t0) ...
                 +  exp( - (x + x0).*(x + x0).*inv4D_t_t0) ...
              ) .*  sqrt(inv4D_t_t0) ./sqrt(pi);  
    



return