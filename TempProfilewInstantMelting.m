% plot some Temperature profile

%Fluence_exp = 2.75;  % kJ/(m^2)
Fluence_exp = 2.2;
Fluence_expsweep = linspace(2.2, 2.9, 4);
for i = 1:7
    
Fluence_exp = Fluence_expsweep(i);
alpha = 2.25e5;     % 1/m
beta = 7e8;         % m fs/kJ

% remove the reflective part
R = 0.33;
Pflux= (1-R)*Fluence_exp;

% critical temperature above which instant (heterogeneous/nonthermal) melting occurs

Ncrit = (5e21)*1e6; % 1/cm^3 (cm^-3 = 1e6 m^-3) 10% electrons
hv = 1.55*1.6e-22; % eV (eV = 1.6e-22 kJ) energy of an electron
Cp = 2410; % kJ/(K m^3) Cp at ~ 1640 K

Tcrit = Ncrit*hv/Cp;
Tmelt = 1685;    %K
%Tcrit = Tmelt;

[TintFct_global,TintInvLam_global] = setupInitTempNonlinAbs_w10electrons(Pflux, beta, alpha);

x = linspace(0,300, 1e4); %in nm the setup function converts to m

ncrit = find(TintFct_global(x)>Tcrit, 1, 'last');

Lv = 4206e3; % crystal latent heat  kJ/m^3

% the resulting temperature in the melt should equilibrate very quickly, so
% it should be the avg temp across the melted region.
Temp1 = quad(TintFct_global,0, x(ncrit))/x(ncrit) - Lv/Cp;

Temp = [Temp1*ones(1, ncrit) TintFct_global(x(ncrit+1:end))];

figure
hold on
% you can turn on and off which temperature profiles you plot by commenting
% Temperature profile if everything thermalized before melting
plot(x, TintFct_global(x), '--k')
% Temperature profile if everything with >10% electrons excited melt (and
% equilibrate across melt) before thermalization
plot(x, Temp)
xlabel('depth [nm]')
ylabel('Temperature [K]')
title(['Fluence   ', sprintf('%0.2g', Fluence_exp),'[kJ/m^2];  nonlinear absorption   ', sprintf('%0.1e',beta), '[m fs/kJ]']) 
savefig(['Fluence', sprintf('%0.2g', Fluence_exp),'nonlinear absorption', sprintf('%0.1e',beta), '.fig'])
end
