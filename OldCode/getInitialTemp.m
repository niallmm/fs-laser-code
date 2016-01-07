function [ TempND, hIC ] = getInitialTemp(Pflux,Tpulse, betain)
global K1 K2
% % Input:
% % Pflux in units kJ/m^2
% % Tpulse in units fs
% % vector of positons in units of 10 nm
% 
% ahatparam = 2e-3; % linear absoption coefficient 
% gammahatparam = 4.5e2;
% c1param = 2e5;
% c2param = 9e7;
% 
% IntND = exp(-ahatparam * xND) ./( 1.0 + PfluxND./TpulseND * gammahatparam * ...
%                                  (1-exp(-ahatparam * xND))                ); 
% 
% DepEnergyND =... % in kJ/m^3
%     PfluxND * c1param * IntND + PfluxND.^2/TpulseND * c2param * IntND .^2;
% 
% 
% MaxDepEnergy = max(DepEnergyND)
% % NOTE THAT SHOULD BE TEMP DEPENDENT
% heatCapND = 2.3*0.88e3; % norm by kJ/m^3/K
% LatHeatND = 4.206e6; % morm. by kJ/m^3
% ThetaND = (1685-300); % in K
% 
% TempSolidND = DepEnergyND./heatCapND./ThetaND;
% TempLiquidND = (DepEnergyND -LatHeatND)./heatCapND./ThetaND;
% %TempLiquidND = max(TempLiquidND,ones(K1+K2,1));
% 

% =================================================
% Niall's version:
%    everything is in m, fs, K, kJ
% =================================================
Tmelt = 1685;    %K
Tambiant = 300;  %K   
Tempin = Tambiant;  % Guess the temp that everything is at
alpha_L = (1.12e5)*exp(Tempin/430);  % 1/m    linear absorption coeff
alpha_FCA = (0.04096)*Tempin/300;      % 1/m??  Free carrier absorption

alpha = alpha_L +alpha_FCA;          % 1/m
%alpha= alphain;

% beta = 9e7; % m fs/kJ


I0 = Pflux/Tpulse; % kJ/(fs m^2)

%for crystaline, and Cp at ~700 K

% Cp = 2024;    % kJ/(K m^3) volumetric heat capacity

Cp = 2410; % Cp at ~ 1640 K

Lv = 4.206e6;    % kJ/m^3 volumetric latent heat
% 
Emelt = Lv + Cp*(Tmelt-Tambiant);
I0Crit = 0.5*(-alpha/beta + sqrt(alpha^2/beta^2 + 4*Emelt/(Tpulse*beta)));



DepEnergy = Tpulse*alpha*alpha*(alpha+beta*I0)*I0*exp(alpha*xND*1e-8)./...
            (beta*I0 - (alpha + beta*I0)*exp(alpha*xND*1e-8)).^2; 
        
% Intense = (alpha/beta)*exp(-alpha*xND*1e-8)./(alpha/(I0*beta) + 1 - exp(-alpha*xND*1e-8));
% 
% DepEnergy = Tpulse*(alpha*Intense+beta*Intense.^2);

        %alpha is in units 1/m and [xND]= 10nm => 10nm*(10e-9 m/nm) = 1e-8 m
        % [DepEnergy] = kJ/m^3
MaxDepEnergy = max(DepEnergy);

        
TempSolidND = DepEnergy./(Cp*(Tmelt-Tambiant)); %no dimensions
TempLiquidND = (DepEnergy-Lv)./(Cp*(Tmelt-Tambiant)); %no dimensions

% MaxTempSolid = max(TempSolidND)*(Tmelt-Tambiant) + Tambiant
% MaxTempLiquid = max(TempLiquidND)*(Tmelt-Tambiant) + Tambiant
MaxTempSolid = max(TempSolidND);
MaxTempLiquid = max(TempLiquidND);

% plot(xND, TempSolidND, '.b')
% hold on
% plot(xND, TempLiquidND, '.r')


%TIC = initialTemp; % ToDothe real factor is the properly rescaled version of (\rho c_p)^{-1}
% we start with a finite melted layer - thus, the temp. in the liquid phase
% is not the one computed above. Thus, one might want to adapt it - in the
% pure equilibrium case it's at melting temp.


TempND(1:K1)   = TempLiquidND(1:K1);
% The temperature should never be smaller than melting temp.... 
% TempND(1:K1) = max(TempND(1:K1),ones(1,K1));

TempND(K1+1) = 0.5 * (TempND(K1) + TempSolidND(K1+1) ); % smoothing
TempND(K1+2:K1+K2) = TempSolidND(K1+2:K1+K2);
TempND(1);
if TempND(1) < 1
    error('Temperature below melting threshold.')
end
if TempND(end) > 1e-9
    error('Temperature profile cut off and energy lost (error tolerance 1e-9), make infinity larger in SetUpParameters.m')
end

