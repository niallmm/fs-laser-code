% parameter sweep
clear all; close all
%Pfluxsweep = linspace(2.5,4,5);
% percentvec = linspace(17, 13, 10);
percentvec = linspace(10, 8.3, 10);
folder = 'percentsweep4';

mkdir(folder);
for sweep = 1:length(percentvec)
    
    % HACKING
    global kappaEglobal vdglobal % only for varying those parameters in the C IRF

    %clear all; close all;
    global xmesh xlabmesh K1 K2 K3 plotflag


    materialpropDir = './checkedCond/';


    addpath(materialpropDir)
    
    SetUpParameters();
    
    %Pflux = 1.9; % in kJ/m^2
%     Pflux = Pfluxsweep(sweep);
Pflux = 1.5;
    tpulse = 100; % in fs
    % we don't know the parameters of the concentration IRF yet - so here one
    % can vary them
    % the function is csolid/cliquid = (ke + hdot/vd )/(1 + hdot/vd)
    % reasonable values for ke are... small (cp. to one) and vd is on the order
    % of severeal m/s (between 1 and 15 m/s - in our units that translates into
    % 0.1 to 1.5

kappaEglobal = 1e-9;
vdglobal     = 0.1;% in 10nm/ns = 10 m/s
 beta = 5.6e8; % m fs/kJ
% beta = 1e9;
percent = percentvec(sweep);
% percent = 13;

    % Data structure u = [T,C,h] = [Tliq,Tsolid,Cliquid,CS,h] (rescaled frame)
%     hIC = 0.1;
    [TIC(1:K1+K2,1), hIC] = getInitialTemp2(Pflux,tpulse,beta,percent, xmesh);
     figure(2)
     semilogx(xmesh*hIC,TIC)
      hold on
    % Note on storing information in the lab frame:
    % We store the actual values on a grid (xlabmesh). For x values larger than
    % Lright = xlabmesh(end) we also set an asymptotic bulk value.
    
    % initial concentration
    
    %concLabIC = 0.0*(sin(0.5 * xlabmesh)).^2; % just for fun!
    concLabIC = 0*xlabmesh;
    % for larger values that Lright the concentration should be a
    % constant
    concLabbulkIC = 0.0;
    
    % initiate vector indicating state of the material
    % 1: liquid, 2: crystalline, 3: amorphous
    stateLabIC = zeros(K3,1);
    stateLabIC = 2 + stateLabIC; % initially everything is a crystalline solid
    stateLabbulkIC = 2; % state for values larger than Lright
    
    tIC = 0.0;
        hICinDriver = hIC
    plotflag = 0;
   [concLabNew,stateLabNew,tend,hverst] = OneShot(hIC,TIC, concLabIC, concLabbulkIC,stateLabIC, stateLabbulkIC,tIC);
 
    filenumber = sprintf('%d', sweep);
    savefile = strcat(folder, '/percent', filenumber, '.mat');
    save(savefile, 'Pflux','concLabNew','stateLabNew','tend','hverst','TIC', 'percent')
    
%     figure(40+sweep)
%     plot(hverst(:,1),hverst(:,2))
%     rmpath(materialpropDir)

 [hmax,hmaxindex] = max(hverst(:,2));
 tmax = hverst(hmaxindex,1)
meltdepth(sweep) = hmax
% 2.3 time to melt
timetomelt(sweep) = tmax
meltingtime_lowerbound(sweep) = tend
resolidvel(sweep) = hmax/(tend-tmax)
hICvec(sweep) = hIC
   
%    close Figure 1
    clear K1 K2 K3 Pflux TIC concLabIC concLabNew concLabbulkIC filenumber ...
        hIC hverst materialpropDir plotflag savefile stateLabIC stateLabNew ...
        stateLabbulkIC sweep tIC tend tpulse xlabmesh xmesh kappaEglobal vdglobal
    
 end
summary = strcat(folder, '/summary_percent.mat');
save(summary, 'percentvec', 'meltdepth', 'timetomelt', ...
    'meltingtime_lowerbound', 'resolidvel', 'hICvec')