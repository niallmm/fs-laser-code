function T0 = findbeta(beta,Pflux)

alpha = 2.25e5;
Tpulse = 100; % fs
I0 = Pflux/Tpulse; % kJ/(fs m^2)
Cp = 2410; % Cp at ~ 1640 K
T0 = Tpulse*(alpha+beta*I0)*I0/(Cp*1385)-1.4;
end