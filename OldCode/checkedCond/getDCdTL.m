function [ DCdTout ] = getDCdTL(Temp, Conc)
    % compute the der. of solute diffusion constant wrt to temp (liquid)
    % IN: Temp(1:K)
    % IN: Conc(1:K)

    K=length(Temp);
    if (length(Conc)~=K)
        error
    end

%    DCdTout = zeros(K,1);

    % Speed that up with a clever vector notation ...
    
%     for ii = 1:K
% 
%         DCdTout(ii,1) = 0.0;
% 
%     end

  %  DCdTout = 5e1+zeros(K,1);%5.e1*Temp;

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
%
% DCdT = -Ea*D0*exp(-Ea/(k*T))/(k*T^2); 
%====================================================

 Tmelt = 1685 ; % melting temperature of silicon
 Tambiant = 300  ; % ambient temperature of the substrate

% D0norm = 3; 
% Ea = 0.1813;  %eV
% k = 8.617e-5; % boltzmann's constant eV/K
% Tmelt = 1685 ; % melting temperature of silicon K
% Tambiant = 300; % ambient temperature of the substrate K
% T = Temp*(Tmelt-Tambiant)+Tambiant;
% dTdTemp = Tmelt-Tambiant;  %term because we are taking 
%                            %derivative with respect to the nondimentional
%                            %temp
%                            
% 
% DCdTout = D0norm*(Ea./(k*(T).^2))*dTdTemp.*...
%     exp(-Ea./(k*(T)));
aparam = 3.0;
bparam = Tambiant * 8.314/17.7e3;
cparam = (Tmelt-Tambiant) * 8.314/17.7e3;

DCdTout = aparam * exp(-1./(bparam + cparam * Temp)) .* cparam ...
    ./(bparam + cparam * Temp)./(bparam + cparam * Temp);



end

