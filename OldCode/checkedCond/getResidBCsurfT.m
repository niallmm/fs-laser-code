function [ResBCsurfT] = getResidBCsurfT(TguessL1,TguessL2,hguess,tglobal);
global dx1
% Boundary conditions at surface (x=0)
% (a) Temperature at the surface - in the documentation \tilde{F} (capital
% F)
% depends on TguessL(2), TguessL(1), hguess, tglobal (explicit time
% dependence 
% Here we need the real data - for now a 'fake' version
%ResBCsurfT = TguessL2-TguessL1; % zero flux
%or
%option (b) constant Temp
%  ResBCsurfT = TguessL1 - 1.0;
Tmelt = 1685 ; % melting temperature of silicon
Tambiant = 300  ; % ambient temperature of the substrate

% CHECK
%ResBCsurfT = 1.0 + hguess*(TguessL2-TguessL1) - TguessL1;

%ResBCsurfT = (TguessL1 - Tambiant)/(Tmelt - Tambiant) - Tambiant; 


% PURE Equlibrium

%(a) fixed heat flux (i.e. insulated boundary - no heat flux
% imposedGradient = 0;% in K/nm 
% L = 1;% nm
% GradientNonD = imposedGradient * L / (Tmelt - Tambiant);
% ResBCsurfT = TguessL2 - TguessL1 - GradientNonD * hguess * dx1;
% % for vanishing gradient (insulated)
% %ResBCsurfT = TguessL2 - TguessL1;

%(b) fixed temperature 
%imposedTemperature = 1900; % in K
%SurfTempNonD = (imposedTemperature - Tambient) / (Tmelt - Tambiant);
%ResBCsurfT = TguessL1 - SurfTempNonD;

% Try imposing the Stefan Bolzman-Law
emmisivity = 0.1;
paramS = emmisivity * 1.5e-9; % should be e-9 - test!!!
delta = Tambiant/(Tmelt - Tambiant);
kappaTherm = getKappaTherm(TguessL1,1); % 1: it's a liquid!!!


ResBCsurfT = -(TguessL2 - TguessL1) + hguess * dx1 * paramS/kappaTherm * ...
    ( TguessL1^4 + 4 * delta * TguessL1^3 + 6 * delta^2 * TguessL1^2 ...
    + 4* delta^3 * TguessL1 + delta^4 );
