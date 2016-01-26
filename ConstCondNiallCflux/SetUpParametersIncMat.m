function [ ] = SetUpParametersIncMat(Tsurface_GTF)
global K1 L1 K2 L2 Z xmesh dx1 dx2 phi didi_1  dwa dwb Aparam Mu1nm
global dt dtmax dhmax convNewton convthreshld maxIterNewton maxSteps deltaTsave min_h_resold
global xlabmesh Lright  K3

% for the Temp Bulk BC BBC
global Dconst_BBC thresh_check_faraway_BBC Tau_BBC deltaTau_BBC TempThreshold_BBC
% for the short time asymptotics STA
global Dconst_STA alphaParam_STA Linfty_STA maxRelError_STA ParamMu_STA  quadgkRelTol_STA


% =========================================================================
% generate mesh
% =========================================================================
% K1 = 400; % no of points in the liquid phase
K1 = 900;
L1 = 1;  % liquid solid boundary in the rescaled coordinates
% K2 = 402; % no of point in the solid phase (for temperature)
  K2 = 904;
L2 = 4e3; % infinity ;-)
%  L2 = 8e3;
Z = K1 + K2;

% [liquid domain, solid domain] - have point x=1 twice to simplify
 xmesh = [linspace(0,L1,K1)';GetSolidGrid(L2,K2,10.5)'];
 dx1 = xmesh(2) - xmesh(1);
 dx2 = xmesh(K1+2:end) - xmesh(K1+1:end-1);
 didi_1 = [NaN;dx2(2:K2-1)./dx2(1:K2-2)];

%check

check = (dx1/dx2(2)>0.8)&&(dx1/dx2(2)<1.2);
%pause

if(~check)
    fprintf('The grid spacing at the moving fornt differs too much between the liquid and the solid region!!!')
    dx1, dx2(1)
    error('error in SetupParametersIncMat')
end
%error

% Weights for one-sided derivative:
dwa=(xmesh(K1+3) + xmesh(K1+2) -2) ./ (xmesh(K1+2) - xmesh(K1+1))./(xmesh(K1+3) - xmesh(K1+1));
dwb=(xmesh(K1+2) + xmesh(K1+1) -2) ./ (xmesh(K1+3) - xmesh(K1+2))./(xmesh(K1+3) - xmesh(K1+1));


% mesh for representing the solute concentration in the solid phase
% (ATTENTION: in lab-coordinates - that means no transformation to x/h(t)
% !!!)
K3 = 1000; % no of datapoints
Lleft = 0; % surface in lab coordinates
Lright = 10 ; % CHANGE!!! % slightly larger than maximal meldting depth

xlabmesh = [linspace(Lleft,Lright,K3)']; % Note: can also be non-uniform - we don't compute derivatives 

% =========================================================================
% Parameters for the integrator
% =========================================================================
dt = 1e-6; % initial timestep
dtmax = 1e-3; % maximal timestep - maybe change that to a maximal delta h
dhmax = 0.01;
phi = 0.51;%0.5+dt; % phi = 1/2: Crank-Nicolson, phi = 1: Backward Euler
convNewton = 1e-10;%1e-10; % convergence threshold for the Newton iteration
convthreshld = 1e-4 ; % for the time stepping check again the 'validity' of the used norm - currently we are using the sup-norm of the residual
maxIterNewton = 2;
maxSteps = 5e5;

% DEBUG
deltaTsave = 5e-4;
% deltaTsave = 2e-2;

min_h_resold = 1e-3; % stop calculation when h is smaller than this value
% NOTE there is an additional hard coded cutpff at 0.99* hIC 

% =========================================================================
% Set up parameters
% =========================================================================
global Tmelt

% Non-dimensional equation accroding to 03-08-2011
% set parameters
%Theta = Tmelt -Tambiant; 
% Tmelt = 1430 for amorphous
 Tmelt = 1685; % for crystal
% Tambiant = 300
% Tmelt - Tambaint = 1130 K for amorphous
% Tmelt - Tambiant = 1385 K for crystal
% kappa_T^\star = 10 J/(s m K)
% \bar{D} = 1e-7;

%L_v from Hoglund paper
% L_V = 2986 amorphous J/cm^3 = 2986e6 J/m^3
% L_V = 4206 crystal   J/cm^3

Aparam = 32.92; % \Theta * \kappa_T^\star / (\bar{D} L_V)
%Aparam = 37.84; % Amorphous


%Tsurface_GTF = 1; % doesn't matter because all diffusion constants are constant
                  % must make actual surface temp from initial condition if
                  % you want to have anything in your Conditions folder
                  % be a function of temp
Mu1nm = 92; % kinetic undercooling coefficient crystaline 
        % changes in getResidIRFmoveT and getJacobIRFmoveT and is needed for 
        % intial temperature Assymptotics in code below
global vdglobal
%         vdglobal = 0.1; % from Aziz fit this is vd = 1 m/s = 0.1 (10s nm)/ns
%         vdglobal = 1; % MUCH more trapping is vd = 10 m/s = 1 (10s nm)/ns

% =========================================================================
% Now everything for the effective temperature BC 
% =========================================================================

% The FIXED Diffusion constant used for x -> infty (physically: Choose
% something for room temperature IF we consider temperature dependent D for
% the full equations)

Tambient = 300;
Dconst_BBC = getDTS(Tambient,2);

% bookkeeping

% for the melting (requires definiting the initial temp and it's
% logarithmic derivative)
% Check if it is evaluated at x=h(t) *  L2 >> sqrt(4 D t)

% set a threshold for h(t)/sqrt(t)

thresh_check_faraway_BBC = 5.0 ....
                   .* 2.*sqrt(Dconst_BBC)./L2;



% for the resolidification part
Tau_BBC = 1e15; % initialize it as super large
deltaTau_BBC = 0.0;

% If the initial temp or the temp at the grid point K1+K2-3
% (resoldification) is smaller than this threshold, the BC is evaluated as
% zero
TempThreshold_BBC = 1e-11;



% =========================================================================
% And all the material parameters for the Temperature Asymptotics
% =========================================================================

% Fixed Diffusion constant for initial time dynamics - what matters is the
% diffusion close to the surface (strong gradients) and we could tjus
% choose a diffusion constant for high temp (> 1 = melting temp)


Dconst_STA = getDTL(Tsurface_GTF,1);
% !!!! Not independent of Aparam
% in non-d form (see 07/27/12)

% reminder: kappatherm measured in units of 10 J/(K m s)

kappatherm = getKappaTherm(Tsurface_GTF,1); % want it in the liquid flag = 1 gives liquid

alphaParam_STA = Dconst_STA ./ kappatherm ./Aparam;

ParamMu_STA = Mu1nm; %make it the same as full code; must be defined right before we call SetUpParam in driver file
% NOTE !!!! cp. to ParamMu in getResidIRFmovT and getJacobIRFmovT


% Numerical parameters
Linfty_STA = 30; % cutoff for the x \to ininity integrals
maxRelError_STA = 1e-8; % Error bound for: abs int_Linfty^{2 Linfty} / int_0^Linfy
quadgkRelTol_STA = 1e-9; % Relative tolerance for the quadgk fct - that was the problem of the 'noise' - the default is only 1e-4 ...
    % took me a LONG time to fin









return

function [xmesh] = GetSolidGrid(xinf,Npoints,param)

 temp = linspace(0,1,Npoints);
 B = param;
 expB = exp(B);
 C = (expB - xinf)/(expB-1);
 A=1-C;
% 
 xmesh = A*exp(B.*temp) +C;

% Play around with another version to alow for more resolution close to the
% enforced BC at xinf

% 
% Nend = 10;
% temp = linspace(0,1,Npoints-Nend);
% B = param;
% expB = exp(B);
% C = (expB - 0.8*xinf)/(expB-1);
% A=1-C;
% 
% xmeshexpPart = A*exp(B.*temp) +C; % last point is at 0.95 xinf
% xmeshend = linspace(0.8*xinf,xinf,Nend+1);
% xmesh = [xmeshexpPart,xmeshend(2:end)];

return
