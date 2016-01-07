function [ ] = SetUpParameters( )
global K1 L1 K2 L2 Z xmesh dx1 dx2 phi didi_1  dwa dwb Aparam
global dt dtmax dhmax convNewton convthreshld maxIterNewton maxSteps deltaTsave
global xlabmesh Lright  K3


% generate mesh
%K1 = 300; % no of points in the liquid phase
K1 = 1000;
L1 = 1;  % liquid solid boundary in the rescaled coordinates
%K2 = 250; % no of point in the solid phase (for temperature)
K2 = 10000;
L2 = 200000; % infinity ;-)
Z = K1 + K2;

% [liquid domain, solid domain] - have point x=1 twice to simplify
%xmesh = [linspace(0,L1,K1)';GetSolidGrid(L2,K2,12)'];
xmesh = [linspace(0,L1,K1)';GetSolidGrid(L2,K2,11)'];
dx1 = xmesh(2) - xmesh(1);
dx2 = xmesh(K1+2:end) - xmesh(K1+1:end-1);
didi_1 = [NaN;dx2(2:K2-1)./dx2(1:K2-2)];

%check
%dx1, dx2(1), dx1/dx2(2)
%error

% Weights for one-sided derivative:
dwa=(xmesh(K1+3) + xmesh(K1+2) -2) ./ (xmesh(K1+2) - xmesh(K1+1))./(xmesh(K1+3) - xmesh(K1+1));
dwb=(xmesh(K1+2) + xmesh(K1+1) -2) ./ (xmesh(K1+3) - xmesh(K1+2))./(xmesh(K1+3) - xmesh(K1+1));


% mesh for representing the solute concentration in the solid phase
% (ATTENTION: in lab-coordinates - that means no transformation to x/h(t)
% !!!)
K3 = 1000; % no of datapoints
Lleft = 0; % surface in lab coordinates
Lright = 60 ; % CHANGE!!! % slightly larger than maximal meldting depth

xlabmesh = [linspace(Lleft,Lright,K3)']; % Note: can also be non-uniform - we don't compute derivatives 

% Parameters for the integrator
dt = 1e-6; % initial timestep
dtmax = 1e0; % maximal timestep - maybe change that to a maximal delta h
dhmax = 0.01;
phi = 0.52;%0.5+dt; % phi = 1/2: Crank-Nicolson, phi = 1: Backward Euler
convNewton = 1e-10;%1e-10; % convergence threshold for the Newton iteration
convthreshld = 1e-4 ; % for the time stepping check again the 'validity' of the used norm - currently we are using the sup-norm of the residual
maxIterNewton = 2;

maxSteps=1e5; % integration steps 

deltaTsave = 2e-5;

% +++++++++++++++++
% Set up parameters
% +++++++++++++++++

% How much graphical output do we want?
%plotflag = 2; % ToDo: implement the output control

% Non-dimensional equation accroding to 03-08-2011
% set parameters
% Aparam = 32.92; % \Theta * \kappa_T^\star / (\bar{D} L_V)

% ==========================================
% Niall's Additional notes on parameters
% ==========================================
Theta = 1685-300; % K
kappa_T = 1e-2; % kJ/(s m K) Thermal Conductivity -- where is this from? 
L_v = 4.206e6;   % kJ/m^3 volumetric latent heat
D_bar = 1e-7;     % diffusion scale = L/tau L = 10nm, tau = 1ns in m^2/s
Aparam = Theta*kappa_T/(D_bar*L_v);
%L_v from Hoglund paper


end

