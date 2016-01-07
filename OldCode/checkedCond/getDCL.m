function [ DCout ] = getDCL(Temp, Conc)
    % compute the solute diffusion constant in the liquid
    % IN: Temp(1:K)
    % IN: Conc(1:K)
%    global tou0 

    K = length(Temp);

    if (length(Conc)~=K)
        error
    end

    %DCout = zeros(K,1);
    
%    DCOut = ones(K,1);
    
%     Dliquid = 3.2e1 ; % sulfur diffusivity in liquid silicon
%     tau0 = 1; % ns
%     L0 = 1; % nm
%     %D0 = tau0/(L0^2);
%     D0 = (L0^2)/tau0;
%     
%     DCout = Dliquid*DCOut./D0;
% 
%     % Speed that up with a clever vector notation ...
    
%     for ii = 1:K
% 
%         DCout(ii,1) = 0.2;
% 
%     end

%DCout = 1.0 + 5e1*Temp + 5e1*Conc.^3;
%DCout = 2.3e-1*ones(K,1);

%====================================================
% Best parameters: DC = D0*exp(-Ea/(k*T)); 
%  D0 = 9.2e-5 m^2/s from Mazur -- probably solid Si
%  D0 = 3e-7 m^2/s  from Tang 2009 liquid
%  Ea = 2.2 eV from Mazur -- probably solid Si
%  Ea = 17.5 kJ/mol == 0.1813 eV from Tang 2009 liquid
%  k = 8.617e-5 eV/K  boltzmann's constant
% 
% Normalize with the following:
%  Temp = (T - Tambiant)/(Tmelt - Tambiant);
%  D0norm = D0/Dbar;
%  Dbar = 1e-7;
%====================================================

 Tmelt = 1685 ; % melting temperature of silicon
 Tambiant = 300  ; % ambient temperature of the substrate

% D0norm = 3; 
% Ea = 0.1813;  %eV
% k = 8.617e-5; % boltzmann's constant eV/K
% Tmelt = 1685 ; % melting temperature of silicon K
% Tambiant = 300  ; % ambient temperature of the substrate K
% T = Temp*(Tmelt-Tambiant)+Tambiant;
% 
% DCout = D0norm*exp(-Ea./(k*(T)));

aparam = 3.0;
bparam = Tambiant * 8.314/17.7e3;
cparam = (Tmelt-Tambiant) * 8.314/17.7e3;

DCout = aparam * exp(-1./(bparam + cparam * Temp));

end

