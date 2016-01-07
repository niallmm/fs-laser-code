function InitTemp  = setupInitTempNonlinAbs_w10electrons(x, Pflux, beta, alpha)

Tmelt = 1685;    %K
% Tmelt = 1430; %K amorphous
Tambiant = 300;  %K
% Tempin = Tambiant;  % Guess the temp that everything is at
% alpha_L = (1.12e5)*exp(Tempin/430);  % 1/m    linear absorption coeff
% alpha_FCA = (0.04096)*Tempin/300;      % 1/m??  Free carrier absorption
%
% alpha = alpha_L +alpha_FCA;          % 1/m
% alpha= alphain;

%beta = 9e7; % m fs/kJ

Tpulse = 100; % fs

I0 = Pflux/Tpulse; % kJ/(fs m^2)

% define critical threshold where 10% electrons are excited.
Ncrit = (5e21)*1e6; % 1/cm^3 (cm^-3 = 1e6 m^-3) 10% electrons
hv = 1.55*1.6e-22; % eV (eV = 1.6e-22 kJ) energy of an electron
Cp = 2410; % kJ/(K m^3) Cp at ~ 1640 K

Tcrit = Ncrit*hv/Cp;

Lv = 4206e3; % crystal latent heat  kJ/m^3

%for crystaline, and Cp at ~700 K

% Cp = 2024;    % kJ/(K m^3) volumetric heat capacity
Cp = 2410; % Cp at ~ 1640 K

x_temp = linspace(0,300, 1e4); %in nm the setup function converts to m

Temp1 = @(x1) (Tpulse*alpha*alpha*(alpha+beta*I0)*I0*exp(alpha*x1*1e-8)./...
    (Cp*(beta*I0 - (alpha + beta*I0)*exp(alpha*x1*1e-8)).^2));

ncrit = find(Temp1(x_temp)>Tcrit, 1, 'last');

hIC = x_temp(ncrit);

Temp2 = quad(Temp1,0, x_temp(ncrit))/x_temp(ncrit) - Lv/Cp;

for i = 1:length(x)
if x(i)<=x_temp(ncrit)
    
    InitTemp(i) = Temp2/(Tmelt - Tambiant);
else
    InitTemp(i) = Temp1(x(i))/(Tmelt-Tambiant);
end
end

    
    
    
%     
%     Temp1 = quad(TintFct_global,0, x(ncrit))/x(ncrit) - Lv/Cp;
%     
%     Temp = [Temp1*ones(1, ncrit) TintFct_global(x(ncrit+1:end))];
%     
%     InitTempFct=(Tpulse*alpha*alpha*(alpha+beta*I0)*I0*exp(alpha*x*1e-8)./...
%         (Cp*(beta*I0 - (alpha + beta*I0)*exp(alpha*x*1e-8)).^2));
%     
    
