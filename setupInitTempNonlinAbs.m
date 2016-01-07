function [InitTempFct, InitTempOneOverLambdaFct] = setupInitTempNonlinAbs(Pflux, beta, alpha)

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

%for crystaline, and Cp at ~700 K

% Cp = 2024;    % kJ/(K m^3) volumetric heat capacity
Cp = 2410; % Cp at ~ 1640 K

InitTempFct=@(x) (Tpulse*alpha*alpha*(alpha+beta*I0)*I0*exp(alpha*x*1e-8)./...
            (Cp*(beta*I0 - (alpha + beta*I0)*exp(alpha*x*1e-8)).^2))/(Tmelt - Tambiant); 
        
InitTempOneOverLambdaFct = @(x)  alpha*(I0*beta + (alpha + I0*beta)*exp(alpha*x*1e-8))*1e-8./ ...
                            ((alpha + beta*I0)*exp(alpha*x*1e-8) - beta*I0);

