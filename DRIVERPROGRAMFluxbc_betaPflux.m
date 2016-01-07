 %Pfluxsweep = linspace(1.921,2.35,5);


%clear all; close all;


foldername = ['test'];%
% Pflux = 2.5*0.6;
% Pfluxsweep = 0.6*2.5*[0.9, 1.1];
% Pfluxsweep = 0.6*2*0.9; %[0.9, 1.1];
mkdir(foldername)
% beta = 4.1355e+08; % interpolated from 50 nm melt depth
%beta = 4.2283e8; % interpolated from 5 ns total melt time
%beta = 4.547e8; % interpolated from 5 ns resolidifcation time
%betasweep = [linspace(5e8, 1e9, 10) linspace(1.5e9, 5e9, 5)];
%Pfluxsweep = 0.6*linspace(1.8, 2.7, 10)
%betasweep = 7.1e8; %, 3e8];

% for sweep2 = 1:length(Pfluxsweep)
   for sweep1 = 1%:length(betasweep)
        beta = 1.3e9;% betasweep(sweep1);
        Pflux = 0.6*2.5;% Pfluxsweep(sweep2)
        sweep2 = 1;

    tic
global xmesh xlabmesh K1 K2 K3 plotflag plotflagIC plotflag_BBC Mu1nm Aparam
plotflag = 0; % Output of the time integratiom: 0: no output, 1: time integration 2: only for debugging
plotflagIC = 0; % Output for the Short Time Asymptotics
plotflag_BBC = 0;% Output for the temperature bulk BC during resolidification


global vdglobal Tmelt cflux

Tmelt = 1685;
% NOTES:
% 1. Idea: have various IRF / BC definitions in different subfolders and
% switch between them by changing the matlab search path - no idea if that
% works or is elegant.....
% materialpropDir = './compareAsymp/';
% materialpropDir = './FullConditions/';
% materialpropDir = './CondConst/';
% materialpropDir = './ConstCondLiqSConcTmelt/';
% materialpropDir = './FullConditions2/';
% materialpropDir = './ConstCondLiq/';
materialpropDir = './ConstCondNiallCFlux/';
addpath(materialpropDir)

vdglobal = 0.1;
cflux = 0;
% cflux = fluxsweep(sweep2);
% setup initial temperature etc
% generate function handles both for the initial temp as a fct of x as well
% as of - d/dx [ log T(x) ]
%alpha = 0; % dummy variable... defined in setupInitTempNonlinAbs.m check this

alpha = 2.25e5;
Tpulse = 100; % fs
I0 = Pflux/Tpulse; % kJ/(fs m^2)
Cp = 2410; % Cp at ~ 1640 K
T0 = Tpulse*(alpha+beta*I0)*I0/(Cp*1385);
L_temp = 1e-9/(alpha+1.9*I0*beta);

global TintFct_global TintInvLam_global
[TintFct_global,TintInvLam_global] = setupInitTempNonlinAbs(Pflux, beta, alpha);

% Mu1nm = 37.0; %amorphous
Mu1nm = 92; % crystaline
% Aparam = 37.84; % Amorphous
Aparam = 32.92; %crystaline
% Aparam = 80;

Tsurface_GTF = TintFct_global(0.0);
SetUpParametersIncMat(Tsurface_GTF);
DT = getDTL(Tsurface_GTF,0);

figure(234)
xtemp = linspace(0,5000,100);
plot(xtemp, log(TintFct_global(xtemp)));
hold on


% Try to 'guess' a good TimeFinal - Idea: to leading order the velocity is
% given by the IRF evaluated at the surface temperature - that gives an
% upper bound for the velocity. The melting threshold should really be
% simply Tsurface >=1

global ParamMu_STA % from the short time asymptotics 
Tsurface_GTF = TintFct_global(0.0);
hdot_GTF = ParamMu_STA .* (Tsurface_GTF - 1.0);
targeth_GTF = 0.02; % This is the target thickness at which e switch - we will get something smaller bc we overestimate the velocity


if (hdot_GTF > 0)
   
    TimeFinal = targeth_GTF./hdot_GTF
    
else
    
    error('below melting threshold')
end

pause(1)


% getTempLabAtFiniteTime returns
% 1. a function handle for a function that computes the Temp. at arbitrary
% points in space 
% 2. the layer thickness h at t=TimeFinal
% 3. Linfty which is the cutoff for the x-> infty integrals in the
% Greens-fct formulation. Note that evaluating the returned function for x
% larger than maybe 0.5 * Linfty becomes problematic
addpath('./shortTimeAsymp/')

[FunctTIC_x_lab,hIC,Linfty] = getTempLabAtFiniteTime(TimeFinal);

% NEW adapt the mesh in lab space
Lright=10 ; % CHANGE!!! % slightly larger than maximal meldting depth

xlabmesh = getXlabMesh(K3,Lright,K1,hIC);


% CHECK parameters.....
% hIC*L2 has to be large enough. Conditions
% A1: >> sqrt(4 D TimeFinal - already coded as threshold for h/sqrt(t)

global thresh_check_faraway_BBC

%hIC./sqrt(TimeFinal) ./ thresh_check_faraway_BBC
check = hIC./sqrt(TimeFinal) <  thresh_check_faraway_BBC;
%pause

if(check)
   hIC./sqrt(TimeFinal) 
   thresh_check_faraway_BBC
   fprintf('Transition time from the asymptitics to the full code TimeFinal is too short. Or L2 too small.') 

   incL2 = thresh_check_faraway_BBC/(hIC./sqrt(TimeFinal) )
   pause(3)
   %   error('Transition time from the asymptitics to the full code TimeFinal is too short. Or L2 too small.') 
end

indexExtrapolate = find(hIC*xmesh <0.5*Linfty,1,'last');
TIC(1:indexExtrapolate,1) = FunctTIC_x_lab(hIC*xmesh(1:indexExtrapolate));
TIC(indexExtrapolate+1:K1+K2) = TintFct_global(hIC*xmesh(indexExtrapolate+1:K1+K2));

% check the initial temperature
% B1: no jump at indexExtrapolate

if(indexExtrapolate < K1+K2)
    
finiteDer = (TIC(indexExtrapolate) - TIC(indexExtrapolate-1)) ./ (xmesh(indexExtrapolate) - xmesh(indexExtrapolate-1) );
jumpDer   = (TIC(indexExtrapolate+1) - TIC(indexExtrapolate)) ./ (xmesh(indexExtrapolate+1) - xmesh(indexExtrapolate) );

check = abs((finiteDer - jumpDer)./finiteDer) < 0.1;

if(~check)
   error('maybe nonsmooth initial temperature') 
end

end
%error('test')

% B2: 'noise in the integrals'
% only observed in the liquid part
maxcurvature = max(abs(diff(diff(TIC(1:K1))))) ;
%meancurvature = mean(abs(diff(diff(TIC(1:K1)))))
mediancurvature = median(abs(diff(diff(TIC(1:K1)))));

check  = (maxcurvature > 1e2 * mediancurvature);
if (check)
    maxcurvature
    mediancurvature
   figure(1)
%   plot(abs(diff(diff(TIC(1:K1)))))
    
    plot(hIC*xmesh(1:K1+3), TIC(1:K1+3));
    title('CHECK: Initial Temp in the liquid - should be smootj!!!');
    xlabel('x');
    ylabel('T');
    
    %   TICsmLiq = smooth(TIC(1:K1),3,'rloess');
%   hold on;
%   plot(abs(diff(diff(TICsmLiq(1:K1)))),'k')
   fprintf('maybe problem with non-smooth initial temperature - check accuracty of the numerical integration in the STA!')
   pause(3)
   %error ('maybe problem with non-smooth initial temperature - check accuracty of the numerical integration in the STA!')
end




% Note on storing information in the lab frame:
% We store the actual values on a grid (xlabmesh). For x values larger than
% Lright = xlabmesh(end) we also set an asymptotic bulk value. 

% initial concentration



IntIC = 1; % amout of surface sulfur (mol / area)
hconcIC = 0.2*hIC;

% Note: This function 
% if sweep >1
%     concLabIC =concLabNew+getSurfaceLayerGaussianApprox(IntIC,hconcIC,xlabmesh);
% else
    concLabIC =getSurfaceLayerGaussianApprox(IntIC,hconcIC,xlabmesh);
% end


% for larger values that Lright the concentration should be a
% constant 
concLabbulkIC = 0.0;

% initiate vector indicating state of the material
% 1: liquid, 2: crystalline, 3: amorphous 
stateLabIC = zeros(K3,1);
stateLabIC = 2 + stateLabIC; % initially everything is a crystalline solid
stateLabbulkIC = 2; % state for values larger than Lright

tIC = TimeFinal;

TimeTest  = 2* TimeFinal ;

[concLabNew,stateLabNew,tend,hverst,ucheck,tcheck, Tsaved, xliquid, Concliqsave] = OneShotTestAdv(hIC,TIC, concLabIC, concLabbulkIC,stateLabIC, stateLabbulkIC,tIC,TimeTest);

% now check the validity of the asymptotics. 
% compare the solution at a time tckeck \approx TimeTest \approx 2 *
% TimeFinal
% with the asymptics for this time
hcheck = ucheck(end);

[FunctTIC_x_labCheck,hIC_asymp,Linfty] = getTempLabAtFiniteTime(tcheck);

% instead of comparing the whole solution, we only compare the melt depth 

asymptRelErr = abs((hcheck - hIC_asymp)./hcheck)
check = asymptRelErr < 1e-2;

if (~check)
   hcheck
   hIC_asymp
   asymptRelErr
   
    figure(98)
    plot(hIC*xmesh(1:K1+K2),TIC(1:K1+K2),'k')
    hold on;
    plot(hcheck*xmesh(1:K1+K2),ucheck(1:K1+K2),'r')
    plot(hIC_asymp*xmesh(1:K1+ceil(K2/2)),FunctTIC_x_labCheck(hIC_asymp*xmesh(1:K1+ceil(K2/2))),'b')
    hold off
    title('CHECK:asymptotics - k at time t, r integrated to 2t, b STA at time 2t, r and b should overlap!!');
    xlabel('x');
    ylabel('T');
   
    fprintf('Relative Error for the STA maybe too large?')
    pause(3)
    %pause
%    error('Relative Error for the STA maybe too large?')
end


%error('test')

% CHECKS AFTER THE CALC:
% 
 figure(44)
 hold on
 plot(hverst(:,1),hverst(:,2))
 title('Meltdepth vs time');
    xlabel('t');
    ylabel('h');
   
 [hmax,hmaxindex] = max(hverst(:,2));
 tmax = hverst(hmaxindex,1)


% 1. Lright slightly larger than hmax
global Lright
check = (hmax <Lright) && (2*hmax > Lright);

if (~check)
   fprintf('check choice of Lright - it should be slightly larger than hmax')
   hmax
   Lright
end

% 2. Info: Compare relevant lenghtsclaes / timescales
% 2.1 minimal depth at which we apply the x-> infty BC
global L2
minDepthBC = hIC*L2
% 2.2 metl depth
meltdepth(sweep2, sweep1) = hmax
% 2.3 time to melt
timetomelt(sweep2, sweep1) = tmax
T0vec(sweep2, sweep1) = T0;
L_tempvec(sweep2, sweep1)= L_temp;

% 2.4 whole meting time -> Yu-Ting's data (NOTE: lower bound bc we don't
% integrate until h=0 -> if accurate number is needed -> extrapolate to h=0
meltingtime_lowerbound(sweep2,sweep1) = tend
resolidvel(sweep2, sweep1) = hmax/(tend-tmax)
%   filenumber2 =sprintf('%d', sweep2);
%   filenumber = sprintf('%d', sweep);
  betan = sprintf('%0.2e', beta);
  fluxn = sprintf('%0.2e', cflux);
  vdn = sprintf('%0.1f',vdglobal);
  Pfluxn = sprintf('%0.2f', Pflux);
  timetomeltn = sprintf('%0.2e', timetomelt);
  IntICn = sprintf('%0.2f', IntIC);
    savefile = strcat(foldername ,'/P',Pfluxn,'beta',betan, '.mat');    
    save(savefile,'Pflux', 'TimeFinal','concLabNew','stateLabNew','tend',...
        'hverst','TIC', 'Tsaved', 'xmesh', 'xlabmesh', 'alpha','beta', ...
         'Mu1nm', 'Aparam', 'DT', 'vdglobal',  'cflux', 'T0', 'L_temp', 'meltdepth')
    
toc
%    end
end
summary = strcat(foldername, '/summary_pflux', fluxn, '.mat');
save(summary,'betasweep', 'Pfluxsweep', 'meltdepth', 'timetomelt', ...
    'meltingtime_lowerbound', 'resolidvel', 'T0vec', 'L_tempvec')
% end
rmpath(materialpropDir)