function [ TempND, hIC ] = getInitialTemp(Pflux,Tpulse,betain,percent, xmesh)
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
beta = betain;

I0 = Pflux/Tpulse; % kJ/(fs m^2)

%for crystaline, and Cp at ~700 K

% Cp = 2024;    % kJ/(K m^3) volumetric heat capacity

Cp = 2410; % Cp at ~ 1640 K

Lv = 4.206e6;    % kJ/m^3 volumetric latent heat

Ncrit = (5e21)*1e6*percent; % 1/cm^3 (cm^-3 = 1e6 m^-3) 10% electrons
hv = 1.55*1.6e-22; % eV (eV = 1.6e-22 kJ) energy of an electron
Cp = 2410; % kJ/(K m^3) Cp at ~ 1640 K

Tcrit = Ncrit*hv/Cp;
% use a fine resolution and temporary xmesh to calculate the locaiton of
% hIC
x_temp = linspace(0,300, 1e4); %in nm the setup function converts to m
% define the Temperature function from absorption
Temp1 = @(x1) (Tpulse*alpha*alpha*(alpha+beta*I0)*I0*exp(alpha*x1*1e-8)./...
    (Cp*(beta*I0 - (alpha + beta*I0)*exp(alpha*x1*1e-8)).^2));
%find location of hIC based on 10% electrons beign excited
ncrit = find(Temp1(x_temp)>Tcrit, 1, 'last');

figure(2) 

semilogx(x_temp, Temp1(x_temp))
hold on
line([1e-1 1e5], [Tcrit Tcrit], 'Color', 'k', 'LineStyle', '--')
% assign hIC value
%hIC = 0.1;%x_temp(ncrit);
hIC = x_temp(ncrit);
% calculate the Temperature in the liquid assuming that we subtract of the
% latent heat of melting and then the temperature equlibrates in the liquid
% very very quickly.
Temp2 = quad(Temp1,0, x_temp(ncrit))/x_temp(ncrit) - Lv/Cp;

% non-dimensionalize the temp in the liquid and the full temp field as if in solid        
TempLiquidND = Temp2/(Tmelt-Tambiant); %no dimensions
TempSolidND = Temp1(xmesh*hIC)/(Tmelt-Tambiant); %no dimensions

figure(3)
semilogx(xmesh*hIC, TempSolidND)
line([1e-1 1e5], [Tcrit Tcrit]/(Tmelt-Tambiant), 'Color', 'k', 'LineStyle', '--')

% assign temp in liquid to liquid mesh points
TempND(1:K1)   = TempLiquidND;

% average last point in liquid with first point in solid for smoothing
TempND(K1+1) = 0.5 * (TempND(K1) + TempSolidND(K1+1) ); % smoothing
TempND(K1+2) = 0.5*(TempND(K1+1) + TempSolidND(K1+2));
TempND(K1+3) = 0.5*(TempND(K1+2) + TempSolidND(K1+3));

% assign temperature field in the solid to correct mesh points.
TempND(K1+3:K1+K2) = TempSolidND(K1+3:K1+K2);
TempND(1);
% figure(4)
% semilogx(xmesh*hIC, TempND)
if TempND(1) < 1
    error('Temperature below melting threshold.')
end
if TempND(end) > 1e-9
    error('Temperature profile cut off and energy lost (error tolerance 1e-9), make infinity larger in SetUpParameters.m')
end

