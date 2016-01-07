function [ DCout ] = getDCL(Temp, Conc)
global Tmelt
% compute the solute diffusion constant in the liquid
    % IN: Temp(1:K)
    % IN: Conc(1:K)


    K = length(Temp);

    if (length(Conc)~=K)
        error
    end



%====================================================
% Best parameters: DC = D0*exp(-Ea/(k*T)); 
%  D0 = 3e-7 m^2/s  from Tang 2009 liquid
%  Ea = 17.5 kJ/mol == 0.1813 eV from Tang 2009 liquid
%  k = 8.617e-5 eV/K  boltzmann's constant
%  D0 = 
% 
% Normalize with the following:
%  Temp = (T - Tambiant)/(Tmelt - Tambiant);
%  D0norm = D0/Dbar;
%  Dbar = 1e-7; m^2/s
%====================================================
% at the melting temperature the diffusion coefficient ~ D0
% then normalize by Dbar to make non-dimensional

DCout = ones(K,1)*2e-1;     % using Aziz fit
% DCout = ones(K,1)*7e-1;
% DCout = ones(K,1)*3;       % using Tang data

end

