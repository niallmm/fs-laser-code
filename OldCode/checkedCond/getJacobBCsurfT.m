function [JBCsurfT] = getJacobBCsurfT(TguessL1, TguessL2, hguess, tglobal)
global K1  K2 Z dx1
% Vector: T[1:K1] liquid, T[K1+1,K1+K2] solid, C[1:K1] liquid, C[K1+1] (SINGLE
% POINT) solid, h

% Boundary conditions at surface (x=0)
% (a) Temperature at the surface - in the documentation \tilde{F} (capital
% F)
% depends on TguessL(2), TguessL(1), hguess, tglobal (explicit time
% dependence 

JBCsurfT = sparse(1,Z+K1+2);
% zero flux

%JBCsurfT(1,2) = 1.0;
%JBCsurfT(1,1) = -1.0;


% %CHECK
% %der. wrt TguessL2
% JBCsurfT(1,2) = hguess;
% % TguessL1
% JBCsurfT(1,1) = -1.0 - hguess;
% 
% % h
% JBCsurfT(1,Z+K1+2) = TguessL2-TguessL1;


Tmelt = 1685 ; % melting temperature of silicon
Tambiant = 300  ; % ambient temperature of the substrate
%JBCsurfT(1,1) = 1/(Tmelt - Tambiant); 

% PURE Equlibrium

%(a) fixed heat flux (i.e. insulated boundary - no heat flux
% imposedGradient = 0;% in K/nm 
% L = 1;% nm
% GradientNonD = imposedGradient * L / (Tmelt - Tambiant);
% %ResBCsurfT = TguessL2 - TguessL1 - GradientNonD * hguess * dx1;
% % TguessL1
% JBCsurfT(1,1) = -1.0;
% % TguessL2
% JBCsurfT(1,2) = 1.0;
% % h
% JBCsurfT(1,Z+K1+2) = -GradientNonD * dx1;

% % for vanishing gradient (insulated)
% %ResBCsurfT = TguessL2 - TguessL1;
% JBCsurgT(1,1) = -1.0;
% % TguessL2
% JBCsurgT(1,2) = 1.0;


%(b) fixed temperature 
%imposedTemperature = 1900; % in K
%SurfTempNonD = (imposedTemperature - Tambient) / (Tmelt - Tambiant);
%ResBCsurfT = TguessL1 - SurfTempNonD;

% TguessL1
%JBCsurgT(1,1) = 1.0;


% Try imposing Stefan-Boltzmann
% NOTE: Check if the radiation is of any importance - the very small
% parameter paramS which is a rescaed version of the stefan constant is
% very small - so maybe an insulatng bc is suffiecient
emmisivity = 0.1;
paramS = emmisivity * 1.5e-9; % should be e-9 - test
delta = Tambiant/(Tmelt - Tambiant);
kappaTherm = getKappaTherm(TguessL1,1); % 1: it's a liquid!!!
dkappaThermdT = getKappaThermdT(TguessL1,1);

% TguessL1
JBCsurfT(1,1) = 1.0 + hguess * dx1 * paramS/kappaTherm * ...
     ( 4 * TguessL1^3 + 12 * delta * TguessL1^2 + 12 * delta^2 * TguessL1 ...
     + 4* delta^3 ...
     ) ... % partial derivative with respect to K_T 
     - dkappaThermdT * ...
        hguess * dx1 * paramS/kappaTherm/kappaTherm * ...
            ( TguessL1^4 + 4 * delta * TguessL1^3 + 6 * delta^2 * TguessL1^2 ...
             + 4* delta^3 * TguessL1 + delta^4 ) ;
% TguessL2
JBCsurfT(1,2) = -1.0;
% h
JBCsurfT(1,Z+K1+2) = dx1 * paramS/kappaTherm * ...
    ( TguessL1^4 + 4 * delta * TguessL1^3 + 6 * delta^2 * TguessL1^2 ...
    + 4* delta^3 * TguessL1 + delta^4 );

        
;
